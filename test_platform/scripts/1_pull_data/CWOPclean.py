### Draft 1: Cleaning script for CWOP network -- solar radiation data ONLY

# import wget
import os
from datetime import datetime, timezone
import gzip
import shutil ## unneccessary now, I think


# Solar radiation data is stored on Google Drive.
# Each folder contains a year of data (2009-2022), and each file
# corresponds to a day of data.

# Step 1: download data

# Step 2: read in and unzip .gz files.
## Can run gunzip from terminal, but find an approach in python first.

#import gzip
#with gzip.open('CWOP/SR/L20180101.txt.gz', 'rb') as f:
#    file_content = f.read()

## Read in .txt file
# date = open("CWOP/SR/L20180101.txt")

# print(date.readline()) # Get line of text to see format.

## Read in CWOP data using gzip

workdir = "/Users/victoriaford/Desktop/era_files/historical_obs/histobs_dev_code/historical-obs-platform/test_platform/scripts/1_pull_data/"

def get_cwop(workdir, get_all = True):
    # for i in ['2021 - 365 files']: #tetsing single year -- google drive folder for SR data
    # public access years: 2018 - 2022
    # private years: 2009 - 2017

    with gzip.open("CWOP/L20210102.txt.gz", "rb") as f_in:
        file_content = f_in.read()
        print(file_content)

        f_in.close() # closes file


        # f_in = open("/CWOP/L20210101.txt")
        # f_out = gzip.open("/CWOP/L20210101.txt.gz", 'wb')
        # f_out.writelines(f_in)
        #
        #
        # # print(date.readline())
        # f_out.close()
        # f_in.close()

        # write code here
        # for filename in filenames:
            # all characters before '>' is station # id number
            # all characters before 'z' is UTC time
            # next 8 characters is latitude followed by 'N'
            # next 9 characters is longitude followed by 'W'
                # if 'E' is present instead, skip lines?
            # '_' 3 characters is wind direction
            # '/' 3 characters plus "g" is wind gusts?
            #

            # if get_all is False: # following maritime file for now
            #     file = open(local_filename, 'wb')
            #     file.close() # closes file
            # elif get_all is True: # if get_all is true, download all new files?? consider google Drive
            #     local_filename = os.path.join(workdir, filename)
            #     file = open(local_filename, 'wb')
            #     file.close()

get_cwop(workdir, get_all = True)
