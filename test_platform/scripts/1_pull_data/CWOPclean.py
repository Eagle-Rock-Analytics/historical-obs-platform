# Draft 1: Cleaning script for CWOP network -- solar radiation data ONLY
# public access years: 2018 - 2022
# private years: 2009 - 2017

import os
from datetime import datetime, timezone
import gzip
import shutil
import zipfile
import xarray as xr

# set envr variables
workdir = "/Users/victoriaford/Desktop/historical_obs/historical-obs-platform/test_platform/scripts/1_pull_data/"
datadir = "/Users/victoriaford/Desktop/historical-obs/historical-obs-platform/test_platform/CWOP_SR/"

# Step 1: download data -- batch download from google drive - folder downloads as a .zip archive
# For new data downloads only need to do once
# Tested on 2020 and 2021 data

# Step 2: unzip .gz batch files -- This works from the command line, but not in script
# unzip drive-download.zip -d datadir   # folder name is dependent on when downloaded
# Step 2b: remove .zip folder from directory?

# def unzip_dir(directory):
#     os.chdir(directory)
#     extension = ".zip"
#     for item in os.listdir(directory):
#         if item.endswith(extension):
#             zip_name = os.path.abspath(item)
#             print('unzipping folder: ', zip_name)
#             unzip zip_name -d directory               # barks at this line
#             os.remove(zip_name)
#
# unzip_dir(datadir)


# Step 3: read in and unzip .gz daily files
def gz_extract(directory): ## modified from https://gist.github.com/kstreepy/a9800804c21367d5a8bde692318a18f5
    extension = '.gz'
    os.chdir(directory)
    for item in os.listdir(directory):
        if item.endswith(extension):
            gz_name = os.path.abspath(item)
            print('unzipping: ', item)
            file_name = (os.path.basename(gz_name)).rsplit('.',1)[0]
            with gzip.open(gz_name, 'rb') as f_in, open(file_name, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(gz_name)  # removes .gz filenames, cleans up directory for ease

gz_extract(datadir)


# Step 4: read in .txt files -- following maritime script
## Code written based off of test file first

#NOTE:  solar data starts at 11pm UTC

def get_cwop(workdir, get_all = True):
    filename = "L20200101.txt"
    date = open(filename)
    print(date.readline())  # Getting line of text to see format

    # file = os.path.join(workdir, filename)
    #
    # # Create a list of variables to remove
    # dropvars = [] # empty for now -- SR data is messy
    # # going to have to cut based off of character limit?
    #
    # # Read in file and drop variables that aren't of interest
    # ds = xr.open_dataset(file, drop_variables=dropvars)
    # print(ds)

get_cwop(datadir, get_all=True)


    # TO DO: configure WD to write files to (in AWS)
