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
from ftplib import FTP
from datetime import datetime, timezone
import pandas as pd
import csv
import boto3
from io import BytesIO, StringIO

# Set AWS credentials
s3 = boto3.client('s3')
bucket_name = 'wecc-historical-wx'
directory = '1_raw_wx/MARITIME/'

# Set envr variables
#workdir = "test_platform/data/1_raw_wx/MARITIME/"
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
    print('{} saved'.format(file)) # Helpful for testing, can be removed.
    r.close() # Close file

# Step 0: Get list of IDs to download.
### Get all stations that start with "46" (Pacific Ocean) and all lettered (all CMAN) stations.

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

                            ### If get_all is False, only download files whose date has changed since the last download.
                            #### Note that the way this is written will NOT fix partial downloads (i.e. it does not check if the specific file is in the folder)
                            #### It will only add new/changed files to a complete set (i.e. add files newer than the last download.)
                            #### This code could be altered to compare write time and file name if desired.
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
    s3.put_object(Bucket=bucket_name, Body=content,Key=directory+"errors_maritime_{}.csv".format(end_api))

    ftp.quit() # This is the “polite” way to close a connection

# To download all data, run:
get_maritime(bucket_name, directory, years = years, get_all = True)

# Note, for first full data pull, set get_all = True
# For all subsequent data pulls/update with newer data, set get_all = False
