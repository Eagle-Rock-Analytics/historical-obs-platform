# Scraping script for MARITIME network.

#### TO DO LIST
# TO DO: specify recipient folder (in AWS?)
## To do so, see code here:
## https://www.linkedin.com/pulse/downloading-files-from-remote-sftp-server-directly-aws-tom-reid (if sftp supported)
## or here: https://stackoverflow.com/questions/64884023/ftps-to-aws-s3-with-ftplib-and-boto3-error-fileobj-must-implement-read (ftplib only)
### Make sure any code with os.dir works when in AWS S3, or adapt.

## Load packages
from ftplib import FTP
import os
from datetime import datetime, timezone
import xarray as xr
import geopandas as gpd
import csv

# Set envr variables
workdir = "test_platform/data/1_raw_wx/MARITIME/"
years = list(map(str,range(1980,datetime.now().year+1))) # Get list of years from 1980 to current year.

# Set up directory to save files, if it doesn't already exist.
try:
    os.mkdir(workdir) # Make the directory to save data in. Except used to pass through code if folder already exists.
except:
    pass


# Step 0: Get list of IDs to download.
### Get all stations that start with "46" (Pacific Ocean) and all lettered (all CMAN) stations.

## Read in MARITIME data using FTP.
def get_maritime(workdir, years, get_all = True):
    
    # ## Login.
    # ## using ftplib
    ftp = FTP('ftp-oceans.ncei.noaa.gov')
    ftp.login() # user anonymous, password anonymous
    ftp.cwd('pub/data.nodc/ndbc/cmanwx/')  # Change WD.
    pwd = ftp.pwd() # Get base file path.

    try:
        os.mkdir(workdir) # Make the directory to save data in. Except used to pass through code if folder already exists.
    except:
        pass
    
    # Set up error handling df.
    errors = {'File':[], 'Time':[], 'Error':[]}

    # Set end time to be current time at beginning of download
    end_api = datetime.now().strftime('%Y%m%d%H%M')

    # Get date of most recently edited file.                   
    try:
        last_edit_time = max([f for f in os.scandir(workdir)], key=lambda x: x.stat().st_mtime).stat().st_mtime
        last_edit_time = datetime.fromtimestamp(last_edit_time, tz=timezone.utc)
    except:
        get_all = True # If folder empty or there's an error with the "last downloaded" metadata, redownload all data.
 
    #for i in years: # For each year/folder
    for i in ['2021']: # Subset for testing (don't download files older than 1980)
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
                    #filenames.append("fake") To test error writing csv.
                    for filename in filenames:
                        try:
                            modifiedTime = ftp.sendcmd('MDTM ' + filename)[4:].strip() # Returns time modified (in UTC)
                            modifiedTime = datetime.strptime(modifiedTime, "%Y%m%d%H%M%S").replace(tzinfo=timezone.utc) # Convert to datetime.
                            
                            ### If get_all is False, only download files whose date has changed since the last download.
                            #### Note that the way this is written will NOT fix partial downloads (i.e. it does not check if the specific file is in the folder)
                            #### It will only add new/changed files to a complete set (i.e. add files newer than the last download.)
                            #### This code could be altered to compare write time and file name if desired.
                            if get_all is False:
                                if (modifiedTime>last_edit_time): # If file new since last run-through, write to folder.
                                    local_filename = os.path.join(workdir, filename) 
                                    file = open(local_filename, 'wb') # Open destination file.
                                    ftp.retrbinary('RETR '+ filename, file.write) # Write file -- EDIT FILE NAMING CONVENTION?
                                    print('{} saved'.format(filename)) # Helpful for testing, can be removed.
                                    file.close() # Close file
                                else:
                                    print("{} already saved".format(filename))
                            elif get_all is True: # If get_all is true, download all files in folder.
                                local_filename = os.path.join(workdir, filename) 
                                file = open(local_filename, 'wb') # Open destination file.
                                ftp.retrbinary('RETR '+ filename, file.write) # Write file -- EDIT FILE NAMING CONVENTION?
                                print('{} saved'.format(filename)) # Helpful for testing, can be removed.
                                file.close() # Close file
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
    filepath = workdir+"errors_cwop_{}.csv".format(end_api) # Set path to save error file.
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


    ftp.quit() # This is the “polite” way to close a connection

# To download data, run:
get_maritime(workdir, years = years, get_all = True)

