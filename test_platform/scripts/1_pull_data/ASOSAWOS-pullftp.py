### Scrape script for ASOS/AWOS network through ISD

# Notes:
# ISD ftp format: Each station has one file per year. The file for each station-year is updated daily for the current year.
# For first pull, use ftp to update station and file.
# To pull real-time data, we may want to write just an API call with date ranges and stations and update the most recent year folder only. 
# This is a separate function/branch.

## Load packages
from ftplib import FTP
import os
from datetime import datetime, timezone
import pandas as pd
from shapely.geometry import shape, Point
import numpy as np
import pandas as pd
import geopandas as gp
from geopandas.tools import sjoin
import boto3 # For AWS integration.
from io import BytesIO

# Set AWS credentials
s3 = boto3.client('s3')
bucket_name = 'wecc-historical-wx'
directory = '1_raw_wx/ASOSAWOS/'

# Set paths to WECC shapefiles in AWS bucket.
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp" 

# Function to return wecc shapefiles and combined bounding box given path variables.
def get_wecc_poly(terrpath, marpath):
    ## get bbox of WECC to use to filter stations against
    ## Read in terrestrial WECC shapefile.
    t = gp.read_file(terrpath)
    ## Read in marine WECC shapefile.
    m = gp.read_file(marpath)
    ## Combine polygons and get bounding box of union.
    bbox = t.union(m).bounds
    return t,m, bbox

# Function to get up to date station list
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
    crs = {'init' :'epsg:4326'} # Set EPSG.
    weccgeo = gp.GeoDataFrame(weccstations, crs=crs, geometry=geometry) # Convert to geodataframe.
    
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
    weccstations['ISD-ID'] = weccstations['USAF']+"-"+weccstations['WBAN'].astype("str")

    # Reformat time strings for FTP/API call.
    weccstations['start_time'] = [datetime.strptime(str(i), '%Y%m%d').strftime('%Y-%m-%d') for i in weccstations['BEGIN']]
    weccstations['end_time'] = [datetime.strptime(str(i), '%Y%m%d').strftime('%Y-%m-%d') for i in weccstations['END']]

    weccstations.reset_index()
    return weccstations

# Function to write FTP data directly to AWS S3 folder.
# ftp here is the current ftp connection
# file is the filename
# directory is the desired path (set of folders) in AWS
def ftp_to_aws(ftp, file, directory):
    r=BytesIO()
    ftp.retrbinary('RETR '+file, r.write)
    r.seek(0)
    s3.upload_fileobj(r, bucket_name, directory+file)
    print('{} saved'.format(file)) # Helpful for testing, can be removed.
    r.close() # Close file

# Function to query ftp server for ISD data. Run this one time to get all historical data or to update changed files for all years.
# Start date format: 'YYYY-MM-DD"
def get_isd_data_ftp(station_list, bucket_name, directory, start_date = None, get_all = True): 
    
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
        s3 = boto3.client('s3')
        objects = s3.list_objects(Bucket=bucket_name,Prefix = directory)
        all = objects['Contents']     
        # Get date of last edited file   
        latest = max(all, key=lambda x: x['LastModified'])
        last_edit_time = latest['LastModified']
        # Get list of all file names
        alreadysaved = []
        for item in all:
            files = item['Key']
            alreadysaved.append(files) 
        alreadysaved = [ele.replace(directory, '') for ele in alreadysaved]
        print(alreadysaved)
    except:
        get_all = True # If folder empty or there's an error with the "last downloaded" metadata, redownload all data.
 
    # Set up AWS to write to bucket.
    s3 = boto3.resource('s3')

    for i in years: # For each year / folder.
    #for i in ['1973', '1989', '2004', '2015', '2021']: # For testing
        if len(i)<5: # If folder is the name of a year (and not metadata file)
            if (start_date is not None and int(i)>int(start_date[0:4])) or start_date is None:  
                # If no start date specified or year of folder is within start date range, download folder.
                try:
                    ftp.cwd(pwd) # Return to original working directory
                    ftp.cwd(i) # Change working directory to year.
                    filenames = ftp.nlst() # Get list of all file names in folder. 
                    filefiltlist = station_list["ISD-ID"]+"-"+i+'.gz' # Reformat station IDs to match file names.
                    filefiltlist = filefiltlist.tolist() # Convert to list.
                    fileswecc = [x for x in filenames if x in filefiltlist] # Only pull all file names that are contained in station_list ID column.
                    #fileswecc = fileswecc[0:5] # For downloading sample of data. FOR TESTING ONLY.
                    for filename in fileswecc:
                        modifiedTime = ftp.sendcmd('MDTM ' + filename)[4:].strip() # Returns time modified (in UTC)
                        modifiedTime = datetime.strptime(modifiedTime, "%Y%m%d%H%M%S").replace(tzinfo=timezone.utc) # Convert to datetime.
                        
                        ### If get_all is False, only download files that are not already in the folder, or that have been updated since the last download.
                        if get_all is False:
                            if filename in alreadysaved: # If filename already in saved bucket
                                if (modifiedTime>last_edit_time): # If file new since last run-through, write to folder.
                                    ftp_to_aws(ftp, filename, directory)    
                                else:
                                    print("{} already saved".format(filename))
                            else:
                                ftp_to_aws(ftp, filename, directory) # Else, if filename not saved already, save.
                            
                        elif get_all is True: # If get_all is true, download all files in folder.
                            ftp_to_aws(ftp, filename, directory) 
                            
                except Exception as e:
                    print("Error in downloading date {}: {}". format(i, e))
                    next  # Adds error handling in case of missing folder. Skip to next folder.
            else: # If year of folder not in start date range, skip folder.
                next
            
        else:
            next # Skip if file or folder isn't a year. Can change to print file/folder name, or to save other metadata files as desired.

    ftp.quit() # This is the “polite” way to close a connection

# Run functions
stations = get_wecc_stations(wecc_terr, wecc_mar)
print(stations)
get_isd_data_ftp(stations, bucket_name, directory, start_date = "1980-01-10", get_all = False)
