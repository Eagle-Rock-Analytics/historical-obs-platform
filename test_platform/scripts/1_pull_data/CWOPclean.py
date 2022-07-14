# Draft 1: Cleaning script for CWOP network -- solar radiation data ONLY
# public access years: 2018 - 2022
# private years: 2009 - 2017
# testing using 2020 and 2021


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
## CWOP data report notes you can connect google drive folder to AWS which may be useful for automation

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
# therefore, each file comprises mix of data from the named day + last hour of previous day
# will need to handle

def get_cwop(workdir, get_all = True):
    filename = "L20200101.txt"
    date = open(filename)
    print(date.readline())  # Getting line of text to see format

    file = os.path.join(workdir, filename)

    # modified from: https://github.com/lohancock/solar-data-parser

    for line in file:  # note this is bad counting on my part for now and just taking rough notes
        line[:5] = 'station_id'
        line[6] = '>'
        line[7:13] = 'time val'
        line[13] = 'z' # utc time
        line[14:21] = 'latitude val' # will need to parse based off of WECC boundary
        line[21] = 'latitude sign'
        line[22] = '/'
        line[23:31] = 'longitude val' # will need to parse based off of WECC boundary
        line[32] = 'longitude sign'
        line[33] = '_'
        line[34:38] = 'wind direction'
        line[38] = '/'
        line[39:42] = 'wind speed val'
        line[42] = 'g' # gust designator
        line[43:46] = 'gust speed val'
        line[47] = 't' # temperature designator
        line[48:51] = 'temperature val'
        line[51] = 'P' # precipitation designator
        line[52:55] = 'precipitation val, units hundreds of inch'
        line[55] = 'h' # relative humidity designator
        line[56:58] = 'relative humidity val, percentage'
        line[58] = 'b' # barometric pressure designator
        line[59:64] = 'air pressure val, units mb with 1 decimal place'
        line[65] = 'L' # irradiance designator
        line[66:69] = 'irradiance val'
        line[69:] = 'tech that was used for rad obs -- undefined length'


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


# TO DO:
# convert date-time to standard format
# convert station metadata to standard format
# lat-lon coordinates
# station elevation -- may need to calculate using R script
# station_id
# network_id
# add dataset metadata in standard format
# data source
# time of download
# citation info
# convert missing data to common format
# convert exisitng qa/qc flags to standard format
