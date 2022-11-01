"""
This script downloads ASOS and AWOS data from ISD using ftp.
Approach:
(1) Get station list (does not need to be re-run constantly)
(2) Download data using station list.
Inputs: bucket name in AWS, directory to save file to (folder path), station list, start date of file pull (optional),
parameter to only download changed files (optional)
Outputs: Raw data for an individual network, all variables, all times. Organized by station, with 1 file per year.

Notes:
1. The file for each station-year is updated daily for the current year. 
To pull real-time data, we may want to write just an API call with date ranges and stations and update the most recent year folder only. 
This is a separate function/branch.
2. This function assumes users have configured the AWS CLI such that their access key / secret key pair are stored in ~/.aws/credentials.
See https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html for guidance.
"""

## Step 0: Environment set-up
# Import libraries
from ftplib import FTP
from datetime import datetime, timezone
import pandas as pd
from shapely.geometry import Point
import pandas as pd
import geopandas as gp
from geopandas.tools import sjoin
import boto3 # For AWS integration.
from io import BytesIO, StringIO
import calc_pull
import requests
import numpy as np

# Set envr variables

# Set AWS credentials
s3 = boto3.client('s3')
bucket_name = 'wecc-historical-wx'
directory = '1_raw_wx/ASOSAWOS/'

# Set paths to WECC shapefiles in AWS bucket.
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp" 

# Set state shortcodes
states = [ 'AK', 'AL', 'AR', 'AZ', 'CA', 'CO', 'CT', 'DC', 'DE', 'FL', 'GA',
           'HI', 'IA', 'ID', 'IL', 'IN', 'KS', 'KY', 'LA', 'MA', 'MD', 'ME',
           'MI', 'MN', 'MO', 'MS', 'MT', 'NC', 'ND', 'NE', 'NH', 'NJ', 'NM',
           'NV', 'NY', 'OH', 'OK', 'OR', 'PA', 'RI', 'SC', 'SD', 'TN', 'TX',
           'UT', 'VA', 'VT', 'WA', 'WI', 'WV', 'WY']

# Function to write FTP data directly to AWS S3 folder.
# ftp here is the current ftp connection
# file is the filename
# directory is the desired path (set of folders) in AWS
def ftp_to_aws(ftp, file, directory):
    r=BytesIO()
    ftp.retrbinary('RETR '+file, r.write)
    r.seek(0)
    s3.upload_fileobj(r, bucket_name, directory+file)
    print('{} saved'.format(file)) # Optional.
    r.close() # Close file

# Function to download and parse ASOS and AWOS station lists (.txt).
# Source: https://www.ncei.noaa.gov/access/homr/reports/platforms
# Inputs: none.
# Outputs: one asosawos-stations.csv file and a stations object.
# This function is called internally in get_wecc_stations().
def get_asosawos_stations():
    
    # Get AWOS stations.
    awosurl = "https://www.ncei.noaa.gov/access/homr/file/awos-stations.txt"
    awosr = requests.get(awosurl)
    lines = awosr.content.split(b'\n') # Split by line
    
    lines = lines[4:] # Skip the first 4 lines
    df = []
    for line in lines: # Parse data manually.
        ncdcid = line[0:8]
        wban = line[9:14]
        coopid = line[15:21]
        call = line[22:26]
        name = line[27:57]
        country = line[58:78]
        st = line[79:81]
        county = line[82:112]
        lat = line[113:122]
        lon = line[123:133]
        elev = line[134:140]
        utc = line[141:146]
        stntype = line[147:197]

        row = [ncdcid, wban, coopid, call, name, country, st, county, lat, lon, elev, utc, stntype]
        row = [x.decode('utf-8') for x in row] # Convert to string
        row = [x.strip() for x in row] # Strip whitespace
        df.append(row)
        
    # # Convert to pandas
    stations = pd.DataFrame(df, columns = ['NCDCID', 'WBAN', 'COOPID', 'CALL', 'NAME', 'COUNTRY', 'ST', 'COUNTY', 'LAT', 'LON', 'ELEV', 'UTC', 'STNTYPE'])
    
    # Get ASOS stations.
    asosurl = "https://www.ncei.noaa.gov/access/homr/file/asos-stations.txt"
    asosr = requests.get(asosurl)
    lines = asosr.content.split(b'\n') # Split by line
    # Skip the first 4 lines
    lines = lines[4:]
    df = []
    for line in lines:
        ncdcid = line[0:8]
        wban = line[9:14]
        coopid = line[15:21]
        call = line[22:26]
        name = line[27:57]
        alt_name = line[58:88]
        country = line[89:109]
        st = line[110:112]
        county = line[113:143]
        lat = line[144:153]
        lon = line[154:164]
        elev = line[165:171]
        utc = line[172:177]
        stntype = line[178:228]
        begdt = line[229:237]
        ghcnd = line[238:249]
        elev_p = line[250:256]
        elev_a = line[257:263]

        row = [ncdcid, wban, coopid, call, name, alt_name, country, st, 
                county, lat, lon, elev, utc, stntype, begdt, ghcnd, elev_p, elev_a]
        row = [x.decode('utf-8') for x in row] # Convert to string
        row = [x.strip() for x in row] # Strip whitespace
        df.append(row)
        
    # # Convert to pandas
    stationsasos = pd.DataFrame(df, columns = ['NCDCID', 'WBAN', 'COOPID', 'CALL', 'NAME', 'ALTNAME', 'COUNTRY', 'ST', 'COUNTY', 'LAT', 'LON', 'ELEV', 'UTC', 'STNTYPE', 'STARTDATE', 'GHCN-DailyID', 'Barometer_elev', 'Anemometer_elev'])
    
    # Now, merge the two dataframes
    asosawosstations = pd.concat([stationsasos, stations], axis = 0, ignore_index=True)

    # Fill any blank spaces with NaN
    asosawosstations = asosawosstations.replace(r'^\s*$', np.nan, regex=True)

    # Drop 6 records without a WBAN (none in WECC)
    asosawosstations = asosawosstations.dropna(subset=['WBAN'])
    # Note: Only 1 record appears in both ASOS and AWOS - Rock Springs AP in WY.
    
    return asosawosstations

# Function to get up to date station list of ASOS AWOS stations in WECC.
# Pulls in ISD station list and ASOSAWOS station list (two separate csvs), joins by WBAN and returns list of station IDs.
# Inputs: path to terrestrial WECC shapefile, path to marine WECC file. 
# Both paths given relative to home directory for git project.
# Outputs: saves 2 sets of metadata to AWS, returns filtered ISD station object for use in get_asosawos_data_ftp().
def get_wecc_stations(terrpath, marpath): #Could alter script to have shapefile as input also, if there's a use for this.
    ## Login.
    ## using ftplib, get list of stations as csv
    filename = 'isd-history.csv'
    ftp = FTP('ftp.ncdc.noaa.gov')
    ftp.login() # user anonymous, password anonymous
    ftp.cwd('pub/data/noaa/')  # Change WD.

    # Read in ISD stations.
    r=BytesIO()
    ftp.retrbinary('RETR '+filename, r.write)
    r.seek(0)
    
    ## Read in csv and only filter to include US stations.
    stations = pd.read_csv(r)
    weccstations = stations[(stations['CTRY']=="US")]

    # Use spatial geometry to only keep points in wecc marine / terrestrial areas.
    geometry = [Point(xy) for xy in zip(weccstations['LON'], weccstations['LAT'])] # Zip lat lon coords.
    weccgeo = gp.GeoDataFrame(weccstations, crs='EPSG:4326', geometry=geometry) # Convert to geodataframe.
    
    ## get bbox of WECC to use to filter stations against
    t, m, bbox = calc_pull.get_wecc_poly(terrpath, marpath) # Call get_wecc_poly.

    # Get terrestrial stations.
    weccgeo = weccgeo.to_crs(t.crs) # Convert to CRS of terrestrial stations.
    terwecc = sjoin(weccgeo.dropna(), t, how='left') # Only keep stations in terrestrial WECC region.
    terwecc = terwecc.dropna() # Drop empty rows.

    # Get marine stations.
    marwecc = sjoin(weccgeo.dropna(), m, how='left') # Only keep stations in marine WECC region.
    marwecc = marwecc.dropna() # Drop empty rows.
    
    # Join and remove duplicates using USAF and WBAN as combined unique identifier.
    weccstations = (pd.concat([terwecc.iloc[:,:11], marwecc.iloc[:,:11]], ignore_index=True, sort =False)
            .drop_duplicates(['USAF','WBAN'], keep='first'))

    # Generate ID from USAF/WBAN combo for API call. This follows the naming convention used by FTP/AWS for file names.
    # Add leading zeros where they are missing from WBAN stations.
    weccstations['ISD-ID'] = weccstations['USAF']+"-"+weccstations['WBAN'].astype("str").str.pad(5, side = "left", fillchar = "0")

    # Reformat time strings for FTP/API call.
    weccstations['start_time'] = [datetime.strptime(str(i), '%Y%m%d').strftime('%Y-%m-%d') for i in weccstations['BEGIN']]
    weccstations['end_time'] = [datetime.strptime(str(i), '%Y%m%d') for i in weccstations['END']]

    # Only keep stations whose end date is after Jan 01, 1980.
    weccstations = weccstations[weccstations['end_time'] >= datetime.strptime("1980-01-01", "%Y-%m-%d")]
    weccstations['end_time'] = [datetime.strptime(str(i), '%Y%m%d').strftime('%Y-%m-%d') for i in weccstations['END']]

    # Now, read in ASOS and AWOS station files and use to filter to only keep ASOS/AWOS stations.
    asosawos = get_asosawos_stations()
    
    # Make columns have the same format
    asosawos['WBAN'] = asosawos['WBAN'].astype(str).str.pad(5,fillchar='0')
    weccstations['WBAN'] = weccstations['WBAN'].astype(str).str.pad(5,fillchar='0')

    # Convert ASOSAWOS elevation to feet
    asosawos['ELEV'] =  asosawos['ELEV'].astype(float) * 0.3048

    # Sort both dataframes by start date
    asosawos = asosawos.sort_values("STARTDATE")
    weccstations = weccstations.sort_values("BEGIN")

    # Only keep ASOS-AWOS stations in WECC.
    # Merging creates duplicates but gives the number of stations expected.
    # For this stage, keep the 2 sets of station metadata separate and just filter using WBAN.
    m1 = weccstations.WBAN.isin(asosawos.WBAN)
    m2 = asosawos.WBAN.isin(weccstations.WBAN)
    
    weccstations = weccstations[m1]
    asosawos = asosawos[m2]

    weccstations.reset_index(inplace = True, drop = True)
    asosawos.reset_index(inplace = True, drop = True)

    # Write ASOS AWOS station list to CSV.
    csv_buffer = StringIO()
    asosawos.to_csv(csv_buffer)
    content = csv_buffer.getvalue()
    s3.put_object(Bucket=bucket_name, Body=content,Key=directory+"stationlist_asosawos.csv")

    # Write filtered ISD station list to CSV.
    csv_buffer = StringIO()
    weccstations.to_csv(csv_buffer)
    content = csv_buffer.getvalue()
    s3.put_object(Bucket=bucket_name, Body=content,Key=directory+"stationlist_isd_asosawos.csv")

    return weccstations

# Function: query ftp server for ASOS-AWOS data and download zipped files.
# Run this one time to get all historical data or to update changed files for all years.
# Inputs: 
# Station_list: Returned from get_wecc_stations() function.
# bucket_name: name of AWS bucket
# directory: folder path within bucket
# Start date: format 'YYYY-MM-DD" (optional)
# get_all: True or False. If False, only download files whose last edit date is newer than
#  the most recent files downloaded in the save folder. Only use to update a complete set of files.
def get_asosawos_data_ftp(station_list, bucket_name, directory, start_date = None, get_all = True): 
    
    # Set up error handling
    errors = {'Date':[], 'Time':[], 'Error':[]}
    end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download

    ## Login.
    ## using ftplib
    ftp = FTP('ftp.ncdc.noaa.gov')
    ftp.login() # user anonymous, password anonymous
    ftp.cwd('pub/data/noaa')  # Change WD.
    pwd = ftp.pwd() # Get base file path.

    # Get list of folders (by year) in main FTP folder.
    years = ftp.nlst()
    
    # Set up AWS to write to bucket.
    s3 = boto3.client('s3')

    # If no start date specified, manually set to be Jan 01 1980.
    if start_date is None:
        start_date = "1980-01-01"
        
    # Remove depracated stations if filtering by time.
    if start_date is not None:
        try:
            station_list = station_list[station_list["end_time"] >= start_date] # Filter to ensure station is not depracated before time period of interest.
        except Exception as e:
            print("Error:", e) # If error occurs here, 
            years = [i for i in years if (len(i)<5 and int(i)>1979)] # function will use years to filter station files.

    try:
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
    except:
        get_all = True # If folder empty or there's an error with the "last downloaded" metadata, redownload all data.

    for i in years: # For each year / folder.
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
                    for filename in fileswecc:
                        modifiedTime = ftp.sendcmd('MDTM ' + filename)[4:].strip() # Returns time modified (in UTC)
                        modifiedTime = datetime.strptime(modifiedTime, "%Y%m%d%H%M%S").replace(tzinfo=timezone.utc) # Convert to datetime.
                        
                        ### If get_all is False, only download files whose last edit date has changed since the last download or whose filename is not in the folder.
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
    csv_buffer = StringIO()
    errors = pd.DataFrame(errors)
    errors.to_csv(csv_buffer)
    content = csv_buffer.getvalue()
    s3.put_object(Bucket=bucket_name, Body=content,Key=directory+"errors_asosawos_{}.csv".format(end_api))
    
if __name__ == "__main__":
    # Run functions
    stations = get_wecc_stations(wecc_terr, wecc_mar)
    get_asosawos_data_ftp(stations, bucket_name, directory, get_all = True)
