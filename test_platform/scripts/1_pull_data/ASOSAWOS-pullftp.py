"""
This script downloads ASOS and AWOS data from ISD using ftp.
Approach:
(1) Get station list (does not need to be re-run constantly)
(2) Download data using station list.
Inputs: path to savedir (directory to save files to), station list (optional), start date of file pull (optional),
parameter to only download changed files (optional)
Outputs: Raw data for an individual network, all variables, all times. Organized by station, with 1 file per year.

Notes:
The file for each station-year is updated daily for the current year. 
To pull real-time data, we may want to write just an API call with date ranges and stations and update the most recent year folder only. 
This is a separate function/branch.
"""

## Step 0: Environment set-up
# Import libraries
from ftplib import FTP
import os
from datetime import datetime, timezone
import pandas as pd
from shapely.geometry import Point
import pandas as pd
import geopandas as gp
import csv
from geopandas.tools import sjoin

# Set envr variables
savedir = "test_platform/data/1_raw_wx/ASOSAWOS/"
wecc_terr = 'test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp'
wecc_mar = 'test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp'    

# Function to return wecc shapefiles and combined bounding box given path variables.
# Inputs: path to terrestrial WECC shapefile, path to marine WECC file. 
# Both paths given relative to home directory for git project.
def get_wecc_poly(terrpath, marpath):
    ## get bbox of WECC to use to filter stations against
    ## Read in terrestrial WECC shapefile.
    t = gp.read_file(terrpath)
    ## Read in marine WECC shapefile.
    m = gp.read_file(marpath)
    ## Combine polygons and get bounding box of union.
    bbox = t.union(m).bounds
    return t,m, bbox

# Function to get up to date station list of ASOS AWOS stations in WECC.
# Pulls in ISD station list and ASOSAWOS station list (two separate csvs), joins by ICAO and returns list of station IDs.
# Inputs: path to terrestrial WECC shapefile, path to marine WECC file. 
# Both paths given relative to home directory for git project.
def get_wecc_stations(terrpath, marpath): #Could alter script to have shapefile as input also, if there's a use for this.
    ## Login.
    ## using ftplib, get list of stations as csv
    filename = 'isd-history.csv'
    ftp = FTP('ftp.ncdc.noaa.gov')
    ftp.login() # user anonymous, password anonymous
    ftp.cwd('pub/data/noaa/')  # Change WD.
    file = open(filename, 'wb') # Open destination file.
    ftp.retrbinary('RETR '+ filename, file.write) # Write file
    print('{} saved'.format(filename))
    file.close()                                 

    ## Read in csv and only filter to include US stations.
    stations = pd.read_csv("isd-history.csv")
    weccstations = stations[(stations['CTRY']=="US")]

    # Use spatial geometry to only keep points in wecc marine / terrestrial areas.
    geometry = [Point(xy) for xy in zip(weccstations['LON'], weccstations['LAT'])] # Zip lat lon coords.
    weccgeo = gp.GeoDataFrame(weccstations, crs='EPSG:4326', geometry=geometry) # Convert to geodataframe.
    
    ## get bbox of WECC to use to filter stations against
    t, m, bbox = get_wecc_poly(terrpath, marpath) # Call get_wecc_poly.

    # Get terrestrial stations.
    weccgeo = weccgeo.to_crs(t.crs) # Convert to CRS of terrestrial stations.
    terwecc = sjoin(weccgeo, t, how='left') # Only keep stations in terrestrial WECC region.
    terwecc = terwecc.dropna() # Drop empty rows.

    # Get marine stations.
    marwecc = sjoin(weccgeo, m, how='left') # Only keep stations in marine WECC region.
    marwecc = marwecc.dropna() # Drop empty rows.
    
    # Join and remove duplicates using USAF and WBAN as combined unique identifier.
    weccstations = (pd.concat([terwecc.iloc[:,:11], marwecc.iloc[:,:11]], ignore_index=True, sort =False)
            .drop_duplicates(['USAF','WBAN'], keep='first'))

    # Generate ID from USAF/WBAN combo for API call. This follows the naming convention used by FTP/AWS for file names.
    # Add leading zeros where they are missing from WBAN stations.
    weccstations['ISD-ID'] = weccstations['USAF']+"-"+weccstations['WBAN'].astype("str").str.pad(5, side = "left", fillchar = "0")

    # Reformat time strings for FTP/API call.
    weccstations['start_time'] = [datetime.strptime(str(i), '%Y%m%d').strftime('%Y-%m-%d') for i in weccstations['BEGIN']]
    weccstations['end_time'] = [datetime.strptime(str(i), '%Y%m%d').strftime('%Y-%m-%d') for i in weccstations['END']]

    # Now, read in ASOSAWOS csv and use to filter to only keep ASOS/AWOS stations.
    # Source: https://www.aviationweather.gov/docs/metar/stations.txt
    # Last downloaded: 08.25.22
    asosawos = pd.read_csv('test_platform/scripts/2_clean_data/asosawos_stations.csv')
    asosawos = asosawos.loc[(asosawos['A']=="A") | (asosawos['A']=="W")] # A = ASOS, W = AWOS
    asosawos['ICAO'] = asosawos['ICAO'].astype(str) # Fix data types
    weccstations['ICAO'] = weccstations['ICAO'].astype(str) # Fix data types
    asosawos = pd.merge(asosawos, weccstations, on = 'ICAO', how = 'inner') # Join by matching ICAO IDs.
    
    asosawos.reset_index()
    return asosawos

# Function: query ftp server for ASOS-AWOS data and download zipped files.
# Run this one time to get all historical data or to update changed files for all years.
# Inputs: 
# Station_list: Returned from get_wecc_stations() function.
# Startdir: path to save directory (relative to top git repository folder)
# Start date: format 'YYYY-MM-DD" (optional)
# get_all: True or False. If False, only download files whose last edit date is newer than
#  the most recent files downloaded in the save folder. Only use to update a complete set of files.
def get_asosawos_data_ftp(station_list, savedir, start_date = None, get_all = True): 
    
    # Set up directory to save files, if it doesn't already exist.
    try:
        os.mkdir(savedir) # Make the directory to save data in. Except used to pass through code if folder already exists.
    except:
        pass

    # Set up error handling
    errors = {'Date':[], 'Time':[], 'Error':[]}
    end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download

    # Remove depracated stations if filtering by time.
    if start_date is not None:
        try:
            station_list = station_list[station_list["end_time"] >= start_date] # Filter to ensure station is not depracated before time period of interest.
        except Exception as e:
            print("Error:", e) # If error occurs here, function will continue without time filtering. Can add "break" to change this to stop code.

    ## Login.
    ## using ftplib
    ftp = FTP('ftp.ncdc.noaa.gov')
    ftp.login() # user anonymous, password anonymous
    ftp.cwd('pub/data/noaa')  # Change WD.
    pwd = ftp.pwd() # Get base file path.

    # Get list of folders (by year) in main FTP folder.
    years = ftp.nlst()
    
    # TO DO: configure WD to write files to (in AWS)
    
    # Get date of most recently edited file. 
    # Note if using AWS may have to change os function to something that can handle remote repositories. Flagging to revisit.
    try:
        last_edit_time = max([f for f in os.scandir(savedir)], key=lambda x: x.stat().st_mtime).stat().st_mtime
        last_edit_time = datetime.fromtimestamp(last_edit_time, tz=timezone.utc)
    except:
        get_all = True # If folder empty or there's an error with the "last downloaded" metadata, redownload all data.
 
    #for i in years: # For each year / folder.
    for i in ['1989', '2004', '2015', '2021']: # For testing
        if len(i)<5: # If folder is the name of a year (and not metadata file)
            if (start_date is not None and int(i)>=int(start_date[0:4])) or start_date is None:  
                # If no start date specified or year of folder is within start date range, download folder.
                try:
                    ftp.cwd(pwd) # Return to original working directory
                    ftp.cwd(i) # Change working directory to year.
                    filenames = ftp.nlst() # Get list of all file names in folder. 
                    filefiltlist = station_list["ISD-ID"]+"-"+i+'.gz' # Reformat station IDs to match file names.
                    filefiltlist = filefiltlist.tolist() # Convert to list.
                    fileswecc = [x for x in filenames if x in filefiltlist] # Only pull all file names that are contained in station_list ID column.
                    fileswecc = fileswecc[0:40] # For downloading sample of data. FOR TESTING ONLY. Comment out otherwise.
                    for filename in fileswecc:
                        modifiedTime = ftp.sendcmd('MDTM ' + filename)[4:].strip() # Returns time modified (in UTC)
                        modifiedTime = datetime.strptime(modifiedTime, "%Y%m%d%H%M%S").replace(tzinfo=timezone.utc) # Convert to datetime.
                        
                        ### If get_all is False, only download files whose date has changed since the last download.
                        #### Note that the way this is written will NOT fix partial downloads (i.e. it does not check if the specific file is in the folder)
                        #### It will only add new/changed files to a complete set (i.e. add files newer than the last download.)
                        #### This code could be altered to compare write time and file name if desired.
                        if get_all is False:
                            if (modifiedTime>last_edit_time): # If file new since last run-through, write to folder.
                                local_filename = os.path.join(savedir, filename) 
                                file = open(local_filename, 'wb') # Open destination file.
                                ftp.retrbinary('RETR '+ filename, file.write) # Write file -- EDIT FILE NAMING CONVENTION?
                                print('{} saved'.format(filename)) # Helpful for testing, can be removed.
                                file.close() # Close file
                            else:
                                print("{} already saved".format(filename))
                        elif get_all is True: # If get_all is true, download all files in folder.
                            local_filename = os.path.join(savedir, filename) 
                            file = open(local_filename, 'wb') # Open destination file.
                            ftp.retrbinary('RETR '+ filename, file.write) # Write file -- EDIT FILE NAMING CONVENTION?
                            print('{} saved'.format(filename)) # Helpful for testing, can be removed.
                            file.close() # Close file
                except Exception as e:
                    print("Error in downloading date {}: {}". format(i, e))
                    errors['Date'].append(i)
                    errors['Time'].append(end_api)
                    errors['Error'].append(e)

                    next  # Adds error handling in case of missing folder. Skip to next folder.
            else: # If year of folder not in start date range, skip folder.
                next
            
        else:
            next # Skip if file or folder isn't a year. Can change to print file/folder name, or to save other metadata files as desired.

    ftp.quit() # This is the “polite” way to close a connection

    #Write errors to csv
    filepath = savedir+"errors_asosawos_{}.csv".format(end_api) # Set path to save error file.
    #print(errors)
    with open(filepath, "w") as outfile:
        # pass the csv file to csv.writer function.
        writer = csv.writer(outfile)

        # pass the dictionary keys to writerow
        # function to frame the columns of the csv file
        writer.writerow(errors.keys())

        # make use of writerows function to append
        # the remaining values to the corresponding
        # columns using zip function.
        writer.writerows(zip(*errors.values()))

# Run functions
stations = get_wecc_stations(wecc_terr, wecc_mar)
print(stations) # For testing.
get_asosawos_data_ftp(stations, savedir, start_date = "1980-01-01", get_all = True)
