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

# Set envr variables
workdir = "/home/ella/Desktop/Eagle-Rock/Historical-Data-Platform/Maritime/"

## Read in MARITIME data using FTP.
def get_maritime(workdir, get_all = True):
    
    ## Login.
    ## using ftplib
    ftp = FTP('ftp-oceans.ncei.noaa.gov')
    ftp.login() # user anonymous, password anonymous
    ftp.cwd('pub/data.nodc/ndbc/cmanwx/')  # Change WD.
    pwd = ftp.pwd() # Get base file path.

    # Get list of folders (by year) in main FTP folder.
    years = ftp.nlst()

    # TO DO: configure WD to write files to (in AWS)
    
    # Get date of most recently edited file.                   
    try:
        last_edit_time = max([f for f in os.scandir(workdir)], key=lambda x: x.stat().st_mtime).stat().st_mtime
        last_edit_time = datetime.fromtimestamp(last_edit_time, tz=timezone.utc)
    except:
        get_all = True # If folder empty or there's an error with the "last downloaded" metadata, redownload all data.
 
    for i in years: # For each year. FOR REAL RUN.
    #for i in ['1973', '1989', '2004', '2015', '2021']: # For testing
        if len(i)<5: # If folder is the name of a year (and not metadata file)
            for j in range(1, 3): # For each month (1-12) ### Return to 13, just 1-3 now for testing.
                try:
                    ftp.cwd(pwd) # Return to original working directory
                    j = str(j).zfill(2) # Pad to correct format
                    dir = ('{}/{}'.format(i,j)) # Specify working directory by year and month
                    #print(dir)
                    ftp.cwd(dir) # Change working directory to year/month.
                    filenames = ftp.nlst() # Get list of all file names in folder. 
                    filenames = filenames[0:5] # For downloading sample of data, DELETE!! for real run!
                    for filename in filenames:
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
                    print("error in downloading date {}/{}: {}". format(j, i, e))
                    next  # Adds error handling in case of missing folder. Skip to next folder.
        else:
            print(i)

    ftp.quit() # This is the “polite” way to close a connection

# To download data, run:
get_maritime(workdir, get_all = False)

