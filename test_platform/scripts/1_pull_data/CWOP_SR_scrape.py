## Scraping script for CWOP network -- solar radiation data only through Google drive
## Notes: Testing using 2020 and 2021
## Notes: Years 2009 thru 2017 will need access
## Notes: Each text file has NO header/metadata -- will be developed in cleaning script

## TO DO: specify data.requests[at]eaglerockanalytics[dot]com as google drive for AWS
## TO DO: Connect google drive folder to AWS which may be useful for automation -- need to test

## Each file is a daily .txt file containing ALL stations globally

import os
import gzip
import shutil

# set envr variables
datadir = "/Users/victoriaford/Desktop/historical-obs/historical-obs-platform/test_platform/CWOP_SR/"

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
