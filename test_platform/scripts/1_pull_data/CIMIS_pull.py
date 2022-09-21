"""
This script downloads CIMIS data from the CADWR using ftp.
Approach:
(1) Get station list (does not need to be re-run constantly)
(2) Download data using station list.
Inputs: bucket name in AWS, directory to save file to (folder path), station list (optional), start date of file pull (optional),
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
import wget

# Set envr variables

# Set AWS credentials
s3 = boto3.client('s3')
bucket_name = 'wecc-historical-wx'
directory = '1_raw_wx/CIMIS/'

# Set paths to WECC shapefiles in AWS bucket.
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp" 

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

# Function to get up to date station list of CIMIS stations.
# Downloads to AWS S3 the Stations List file most recently found in FTP server.
# Inputs: file path in AWS bucket.
# Note all stations in CA, so no need to filter out stations here. File just downloaded for reference.
def get_cimis_stations(directory): #Could alter script to have shapefile as input also, if there's a use for this.
    ## Login.
    ## using ftplib, get list of stations as csv
    filename = 'CIMIS Stations List (January20).xlsx'
    ftp = FTP('ftpcimis.water.ca.gov')
    ftp.login() # user anonymous, password anonymous

    ftp.cwd('pub2')  # Change WD.
    ftp_to_aws(ftp, filename, directory)
    return
    
# Function: query ftp server for CIMIS data and download csv files.
# Run this one time to get all historical data or to update changed files for all years.
# Inputs: 
# bucket_name: name of AWS bucket
# directory: folder path within bucket
# years: optional, specify years to download data from as a range of years mapped to strings (e.g. ['1998', '1999'])
# get_all: True or False. If False, only download files whose last edit date is newer than
#  the most recent files downloaded in the save folder. Only use to update a complete set of files.
def get_cimis_data_ftp(bucket_name, directory, years = None, get_all = True): 
    
    # Set up error handling
    errors = {'File':[], 'Time':[], 'Error':[]}
    end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download

    ## Login.
    ## using ftplib
    ftp = FTP('ftpcimis.water.ca.gov')
    ftp.login() # user anonymous, password anonymous
    ftp.cwd('pub2/annual')  # Change WD.

    # Set up AWS to write to bucket.
    s3 = boto3.client('s3')

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
        #print(alreadysaved)
    except:
        get_all = True # If folder empty or there's an error with the "last downloaded" metadata, redownload all data.

    try:
        filenames = ftp.nlst() # Get list of all file names in folder. 
        filenames = [i for i in filenames if i.endswith('.zip')]
        filenames = [i for i in filenames if i.startswith('hourly')] # Only keep hourly data, drop daily files.
        if years is not None:
            filenames = [str for str in filenames if any(sub in str for sub in years)] # Filter files by year.
        for filename in filenames:
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
        print("Error in downloading file {}: {}". format(filename, e))
        errors['File'].append(filename)
        errors['Time'].append(end_api)
        errors['Error'].append(e)   


    ftp.quit() # This is the “polite” way to close a connection

    #Write errors to csv
    csv_buffer = StringIO()
    errors = pd.DataFrame(errors)
    errors.to_csv(csv_buffer)
    content = csv_buffer.getvalue()
    s3.put_object(Bucket=bucket_name, Body=content,Key=directory+"errors_cimip_{}.csv".format(end_api))
    
# Run functions
get_cimis_stations(directory)
get_cimis_data_ftp(bucket_name, directory, get_all = True)
# Tested with years = ['1990'] -> works
# Tested with get_all = False -> works
