## Scraping script for CWOP network -- solar radiation data only through Google drive
## Notes: Testing using 2020 and 2021
## Notes: Years 2009 thru 2017 will need access
## Notes: Each text file has NO header/metadata -- will be developed in cleaning script

## TO DO: specify data.requests[at]eaglerockanalytics[dot]com as google drive for AWS
## TO DO: Connect google drive folder to AWS which may be useful for automation -- need to test

## Each file is a daily .txt file containing ALL stations globally

import os
import gzip
from zipfile import ZipFile
import shutil
import re

# set envr variables
homedir = os.getcwd() # Get current working directory.
if "historical-obs-platform" in homedir: # If git folder in path
    homedir = homedir[0:homedir.index("historical-obs-platform")]+"historical-obs-platform" # Set path to top folder.
    os.chdir(homedir) # Change directory.
else:
    print("Error: Set current working directory to the git repository or a subfolder, and then rerun script.")
    exit()

raw_datadir = homedir + "/test_platform/data/1_raw_wx/CWOP_SR/"

## Set up directory to save files, if it doesn't already exist.
# try:
#     os.mkdir(raw_datadir) # Make folder to save cleaned data
#     print("Directory for {} created".format(re.split("/", raw_datadir)[-2]))
# except:
#     print("Directory for {} exists".format(re.split("/", raw_datadir)[-2]))
#     pass    # Pass if folder already exists


## Step 1: download data -- batch download from google drive - folder downloads as a .zip archive
## For new data downloads only need to do this once

## Step 2a: unzip .gz batch files -- This works from the command line, but not in script
# unzip drive-download-20220713T204303Z-001.zip -d datadir   # zip folder name is dependent on when downloaded

## Step 2b: remove .zip folder from directory -- IN PROGRESS
def unzip_dir(directory):
    os.chdir(directory)
    extension = ".zip"
    for item in os.listdir(directory):
        if item.endswith(extension):
            zip_name = os.path.abspath(item)
            print('Unzipping year (folder): ', zip_name)
            ZipFile(zip_name).extractall(directory) # This unzips all files to the directory. (Note you should download the files, rather than the folders from GDrive as the zip will include the file structure.)
            # unzip zip_name -d directory               # barks at this line in script
            os.remove(zip_name)

unzip_dir(datadir)

## Step 3: read in and unzip .gz daily files
## Modified from https://gist.github.com/kstreepy/a9800804c21367d5a8bde692318a18f5
def gz_extract(directory):
    extension = '.gz'
    os.chdir(directory)
    for item in os.listdir(directory):
        if item.endswith(extension):
            gz_name = os.path.abspath(item)
            file_name = (os.path.basename(gz_name)).rsplit('.',1)[0]
            with gzip.open(gz_name, 'rb') as f_in, open(file_name, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(gz_name)  # Removes .gz filenames, cleans up directory for ease
            print('{} saved'. format(file_name)) # Useful check, but can be cleaned up

gz_extract(datadir)

# Useful to add: Error messages/already downloaded messages?
