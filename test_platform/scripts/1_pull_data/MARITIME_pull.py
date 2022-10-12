"""
This script downloads MARITIME data from NCEI using ftp.
Approach:
(1) Download data using station list.
Inputs: bucket name in AWS, directory to save file to (folder path), range of years of data to download,
parameter to only download changed files (optional)
Outputs: Raw data for an individual network, all variables, all times. Organized by station, with 1 file per month.
Notes:
1. This function assumes users have configured the AWS CLI such that their access key / secret key pair are stored in ~/.aws/credentials.
See https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html for guidance.
"""
## Load packages
import requests
from ftplib import FTP
from datetime import datetime, timezone
import pandas as pd
import csv
import boto3
from io import BytesIO, StringIO
import calc_pull
from shapely.geometry import Point
import geopandas as gp
from geopandas.tools import sjoin


# Set AWS credentials
s3 = boto3.client('s3')
s3_cl = boto3.client('s3') # for lower-level processes
bucket_name = 'wecc-historical-wx'
directory_mar = '1_raw_wx/MARITIME/'
directory_ndbc = '1_raw_wx/NDBC/'

# Set paths to WECC shapefiles in AWS bucket.
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"

# Set envr variables
years = list(map(str,range(1980,datetime.now().year+1))) # Get list of years from 1980 to current year.

# Function to write FTP data directly to AWS S3 folder.
# ftp here is the current ftp connection
# file is the filename
# directory is the desired path (set of folders) in AWS
def ftp_to_aws(ftp, file, directory):
    r=BytesIO()
    ftp.retrbinary('RETR '+file, r.write)
    r.seek(0)
    s3.upload_fileobj(r, bucket_name, directory+file)
    print('{}: {} saved'.format(directory, file)) # Helpful for testing which directory data is being saved to, can be removed
    r.close() # Close file

# Step 0: Get list of IDs to download.
### Get all stations that start with "46" (Pacific Ocean) and all lettered (all CMAN) stations.
def get_maritime_station_ids(terrpath, marpath):
    """Returns list of station ids as CSV"""
    url = 'https://www.ndbc.noaa.gov/data/stations/station_table.txt'
    r = requests.get(url)
    lines = r.content.split(b'\n') # split by line
    df = []
    for line in lines:
        row = line.split(b'|') # Split on |
        row = [x.decode('utf-8') for x in row] # converts to string
        row = list(map(str.strip, row))
        df.append(row)

    stations = pd.DataFrame(columns=['STATION_ID', 'OWNER', 'TTYPE', 'HULL', 'NAME', 'PAYLOAD', 'LOCATION', 'TIMEZONE', 'FORECAST', 'NOTE'], data=df[1:])
    stations = stations.dropna()
    stations = stations.drop(columns=['TTYPE', 'HULL', 'PAYLOAD', 'FORECAST', 'TIMEZONE'])

    ## Use spatial geometry to only keep points within WECC marine/terrestrial region (should only be marine)
    ## Locations listed as one value: "30.000 N 90.000 W" as an example, split on spaces
    mar_lon = []
    mar_lat = []
    location = stations['LOCATION']
    for line in location:
        line = line.split(' ')
        if "S" in line:
            mar_lat.append(-1 * float(line[0]))
        else:
            mar_lat.append(float(line[0]))

        if "E" in line:
            mar_lon.append(float(line[2]))
        else:
            mar_lon.append(-1 * float(line[2]))

    stations['LATITUDE'] = mar_lat
    stations['LONGITUDE'] = mar_lon

    geometry = [Point(xy) for xy in zip(stations['LONGITUDE'], stations['LATITUDE'])] # Zip lat lon coords
    weccgeo = gp.GeoDataFrame(stations, crs="EPSG:4326", geometry=geometry) # converts to geodataframe

    ## get bbox of WECC to use to filter stations against
    t, m, bbox = calc_pull.get_wecc_poly(terrpath, marpath) # Call get_wecc_poly.

    # Get terrestrial stations.
    weccgeo = weccgeo.to_crs(t.crs) # Convert to CRS of terrestrial stations.
    terwecc = sjoin(weccgeo.dropna(), t, how='left') # Only keep stations in terrestrial WECC region.
    terwecc = terwecc.dropna() # Drop empty rows.

    # Get marine stations.
    marwecc = sjoin(weccgeo.dropna(), m, how='left') # Only keep stations in marine WECC region.
    marwecc = marwecc.dropna() # Drop empty rows.

    weccstations = (pd.concat([terwecc, marwecc], ignore_index = True, sort=False)
        .drop_duplicates(['STATION_ID'], keep='first'))

    # drop columns and rename in_wecc to in_terr
    weccstations.drop(['OBJECTID_1', 'OBJECTID', 'Shape_Leng','geometry', 'FID_WECC_B', 'BUFF_DIST', 'index_right'], axis = 1, inplace = True)
    weccstations.rename(columns = {'in_WECC':'in_terr_wecc', 'in_marine':'in_mar_wecc'}, inplace = True)

    # Identify which buoys are NDBC and which are CMAN/MARITIME
    network = []
    for item in weccstations['STATION_ID']:
        if item[:2] == "46" or item[:4] == "NDBC_46":
            network.append('NDBC')
        else:
            network.append('MARITIME')

    weccstations['NETWORK'] = network

    # Splits dataframe into respective networks
    maritime_network = weccstations[weccstations['NETWORK'] == 'MARITIME']
    ndbc_network = weccstations[weccstations['NETWORK'] == 'NDBC']

    ## Write stations to respective AWS bucket - requires separate buffers
    mar_buffer = StringIO()
    maritime_network.to_csv(mar_buffer)
    content = mar_buffer.getvalue()
    s3_cl.put_object(Bucket=bucket_name, Body=content, Key=directory_mar+"MARITIME_stations.csv")

    ndbc_buffer = StringIO()
    ndbc_network.to_csv(ndbc_buffer)
    content = ndbc_buffer.getvalue()
    s3_cl.put_object(Bucket=bucket_name, Body=content, Key=directory_ndbc+"NDBC_stations.csv")
    return stations

get_maritime_station_ids(wecc_terr, wecc_mar)


## Read in MARITIME data using FTP.
def get_maritime(bucket_name, directory, years, get_all = True):

    # ## Login.
    # ## using ftplib
    ftp = FTP('ftp-oceans.ncei.noaa.gov')
    ftp.login() # user anonymous, password anonymous
    ftp.cwd('pub/data.nodc/ndbc/cmanwx/')  # Change WD.
    pwd = ftp.pwd() # Get base file path.

    # Set up error handling df.
    errors = {'File':[], 'Time':[], 'Error':[]}

    # Set end time to be current time at beginning of download
    end_api = datetime.now().strftime('%Y%m%d%H%M')

    # Get date of most recently edited file and list of file names already saved.
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

    except:
        get_all = True # If folder empty or there's an error with the "last downloaded" metadata, redownload all data.

    for i in years: # For each year/folder
    #for i in ['2021']: # Subset for testing (don't download files older than 1980)
        if len(i)<5: # If folder is the name of a year (and not metadata file)
            for j in range(1, 13): # For each month (1-12)
                try:
                    ftp.cwd(pwd) # Return to original working directory
                    j = str(j).zfill(2) # Pad to correct format
                    dir = ('{}/{}'.format(i,j)) # Specify working directory by year and month
                    #print(dir)
                    ftp.cwd(dir) # Change working directory to year/month.
                    filenames = ftp.nlst() # Get list of all file names in folder.
                    filenames = list(filter(lambda f: f.endswith('.nc'), filenames)) # Only keep .nc files.
                    filenames = list(filter(lambda f: (f.startswith('46') or f.startswith('NDBC_46') or (f[0].isalpha() and not f.startswith("NDBC"))), filenames)) # Only keep files with start with "46" (Pacific Ocean) or a letter (CMAN buoys.)
                    #filenames.append("fake") #To test error writing csv.
                    for filename in filenames:
                        try:
                            modifiedTime = ftp.sendcmd('MDTM ' + filename)[4:].strip() # Returns time modified (in UTC)
                            modifiedTime = datetime.strptime(modifiedTime, "%Y%m%d%H%M%S").replace(tzinfo=timezone.utc) # Convert to datetime.

                            ## Sorts stations based on whether it is a NDBC buoy or in MARITIME/CMAN
                            if filename[:2] == "46" or filename[:7] == "NDBC_46":
                                if get_all is False:
                                    if filename in alreadysaved: # If filename already in saved bucket
                                        if (modifiedTime>last_edit_time): # If file new since last run-through, write to folder.
                                            ftp_to_aws(ftp, filename, directory=directory_ndbc)
                                        else:
                                            print("{} already saved".format(filename))
                                    else:
                                        ftp_to_aws(ftp, filename, directory=directory_ndbc) # Else, if filename not saved already, save.
                                elif get_all is True: # If get_all is true, download all files in folder.
                                    ftp_to_aws(ftp, filename, directory=directory_ndbc)
                            else:
                                if get_all is False:
                                    if filename in alreadysaved: # If filename already in saved bucket
                                        if (modifiedTime>last_edit_time): # If file new since last run-through, write to folder.
                                            ftp_to_aws(ftp, filename, directory=directory_mar)
                                        else:
                                            print("{} already saved".format(filename))
                                    else:
                                        ftp_to_aws(ftp, filename, directory=directory_mar) # Else, if filename not saved already, save.
                                elif get_all is True: # If get_all is true, download all files in folder.
                                    ftp_to_aws(ftp, filename, directory=directory_mar)
                            ### If get_all is False, only download files whose date has changed since the last download.
                            #### Note that the way this is written will NOT fix partial downloads (i.e. it does not check if the specific file is in the folder)
                            #### It will only add new/changed files to a complete set (i.e. add files newer than the last download.)
                            #### This code could be altered to compare write time and file name if desired.

                        except Exception as e:
                            errors['File'].append(filename)
                            errors['Time'].append(end_api)
                            errors['Error'].append(e)
                except Exception as e:
                    print("error in downloading date {}/{}: {}". format(j, i, e))
                    next  # Adds error handling in case of missing folder. Skip to next folder.
        else:
            print(i)

    # Write errors to csv
    csv_buffer = StringIO()
    errors = pd.DataFrame(errors)
    errors.to_csv(csv_buffer)
    content = csv_buffer.getvalue()
    s3.put_object(Bucket=bucket_name, Body=content,Key=directory_mar+"errors_maritime_{}.csv".format(end_api))

    ftp.quit() # This is the “polite” way to close a connection

# To download data, run:
get_maritime(bucket_name, directory, years = years, get_all = True)
