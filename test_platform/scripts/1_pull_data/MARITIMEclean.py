# Draft 1: Cleaning script for MARITIME network.

#### TO DO LIST
# TO DO: specify recipient folder (in AWS?)
## To do so, see code here:
## https://www.linkedin.com/pulse/downloading-files-from-remote-sftp-server-directly-aws-tom-reid (if sftp supported)
## or here: https://stackoverflow.com/questions/64884023/ftps-to-aws-s3-with-ftplib-and-boto3-error-fileobj-must-implement-read (ftplib only)

# To do: when downloading, compare date of last download to date of file change to only download changed files.
# write or store download date somehow
# Call it in order to compare to new file changes.


## Load packages
from ftplib import FTP
from pydap.client import open_url
import wget
import os
from datetime import datetime, timezone

## Read in MARITIME data using FTP.


## Login.
## using ftplib
ftp = FTP('ftp-oceans.ncei.noaa.gov')
ftp.login() # user anonymous, password anonymous
ftp.cwd('pub/data.nodc/ndbc/cmanwx/')  # Change WD.
pwd = ftp.pwd() # Get base file path.

# Get list of folders (by year) in main FTP folder.
years = ftp.nlst()

# TO DO: configure WD to write files to (in AWS)

# TO FIX: Get date of most recently edited file.
last_edit_time = max([f for f in os.scandir("/home/ella/Desktop/Eagle-Rock/Historical-Data-Platform/Maritime/")], key=lambda x: x.stat().st_mtime).stat().st_mtime
last_edit_time = datetime.fromtimestamp(last_edit_time, tz=timezone.utc)
print(last_edit_time)


#for i in years: # For each year. FOR REAL RUN.
for i in ['1973', '1989', '2004', '2015', '2021']: # For testing
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
                    modifiedTime = datetime. ### FINISH CONVERSION HERE.
                    print(filename)
                    print(modifiedTime)
                    # To do: isolate YD and compare to 
                    ### add toggle in function above, if grab_all = FALSE then use this filtering function.
                    if (modifiedTime>last_edit_time): # If file new since last run-through, write to folder.
                        local_filename = os.path.join('/home/ella/Desktop/Eagle-Rock/Historical-Data-Platform/Maritime/', filename) # Update this with desired data storage location (AWS server)
                        file = open(local_filename, 'wb') # Open destination file.
                        ftp.retrbinary('RETR '+ filename, file.write) # Write file -- EDIT FILE NAMING CONVENTION?
                        print('{} saved'.format(filename))
                        file.close() # Close file

                # For each file, open and drop variables that aren't of interest here, then save again?

            except Exception as e:
                print("error in downloading date {}/{}: {}". format(j, i, e))
                next  # Add error handling in case of missing folder. Skip to next folder.
    else:
        print(i)



ftp.quit() # This is the “polite” way to close a connection
