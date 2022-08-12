"""Draft cleaning script for CWOP network -- solar radiation data ONLY"""

"""
This script is a template structure for data cleaning for a variety of data sources for
ingestion into the Historical Observations Platform.
Approach:
(1) Read through variables, and calculates derived priority variables if not observed
(2) Drops unnecessary variables
(3) Converts station metadata to standard format, with unique identifier
(4) Converts data metadata to standard format, and converts units into standard units if not provided in standard units.
(5) Converts missing data to standard format
(6) Tracks existing qa/qc flag for review
(7) Merge files by station, and outputs cleaned variables as a single .nc file for an individual network.
Inputs: Raw data for an individual network
Outputs: Cleaned data for an individual network, priority variables, all times. Organized by station as .nc file.
"""

# Testing using 2020 and 2021

## Step 0: Environment set-up
## Import libraries
import os
from datetime import datetime, timezone
import xarray as xr
import re
import pandas as pd
# from .thermo_helpers import f_to_c

# set envr variables

datadir = "../../../../historical-obs-platform/test_platform/CWOP_SR/"

## Step 1: Read through variables, and calculates variables if not observed

# Hardcoding boundaries in for now (can use WECC poly in future)
lonmin, lonmax = -139.047795, -102.03721
latmin, latmax = 30.142739, 60.003861

filename = "L20200101.txt"
file = os.path.join(datadir, filename) # pathway to file
date = open(file)
print(date.readline())  # Getting line of text to see format


def parse_cwop_sr(filepath):
    """
    Parses solar radiation data for CWOP.
    Paramters: filepath (str): filepath for file to be parsed
    Returns: data (xr.array): parsed raw data
    """
    lines = []
    date = open(filepath)
    f = open(filepath, 'r')
    count = 0

    while True:
        count += 1
        line = f.readline()

        # Strips off the leading station_id and the trailing hardware specifications, both of which are not of equal character length
        line_items_step1 = re.split(r"[>]", line) # Strip off station_id first
        station_id = line_items_step1[0]
        hardware = line_items_step1[1][66:]     # Hardware specs saved as string -- avoids having to maintain list of all available hardware options (which are many)

        line_items_step2 = re.split(r"[z/_gt]", line_items_step1[1][:66])   # Strips part of remaining string; note: after "t" the variable order can vary per station

        ## Removing data outside WECC region first before any other conversions
        # Lat-lon conversion: stip hemisphere designator, convert to decimal degrees, removes data outside WECC region of interest
        if line_items_step2[1][7] == "N":
            lat_raw = float(line_items_step2[1][:-1]) # Strips the N char off and converts to float
            deg = float(line_items_step2[1][:2])
            min_ = float(line_items_step2[1][2:4])
            sec_ = float(line_items_step2[1][5:7])
            lat_clean = deg + min_/60 + sec_/3600
            if lat_clean < latmin: # getting rid of latitude locations outside WECC region (hard coded)
                continue
            elif lat_clean > latmax:
                continue
            else:
                lat_clean = lat_clean
        elif line_items_step2[1][7] == "S":
            continue
        else:   # primarily error handling in case a line does not have N or S provided -- is this the case or can we get rid of?
            continue
            print("Error: this latitude does not exist/incorrect or missing hemisphere designator.")

        if line_items_step2[2][8] == "W":
            lon_raw = float(line_items_step2[2][:-1]) # Strips the W char off and converts to float
            deg = float(line_items_step2[2][:3])
            min_ = float(line_items_step2[2][3:5])
            sec_ = float(line_items_step2[2][6:8])
            lon_clean = (deg + min_/60 + sec_/3600) * -1 # converting to negative values for standard western hemi reporting
            if lon_clean < lonmin:
                continue
            elif lon_clean > lonmax:
                continue
            else:
                lon_clean = lon_clean
        elif line_items_step2[2][8] == "E":
            continue
        else:   # primarily error handling in case a line does not have E/W provided -- is this the case or can we get rid of?
            continue
            print("Error: this longitude does not exist/incorrect or missing hemisphere designator.")


### give this some more thought -- specifically this is in utc time
### This is working for the first day of the month, but not otherwise

        ###### Note: Solar data starts at 11pm UTC in each file - each file comprises mix of data from the named day + last hour of previous day
        ## Each daily file starts at 11pm UTC the day prior -- will need flag this properly
        ## There are also bad data in some files
        filename_date_raw = datetime.strptime(filename[1:-4].strip(), "%Y%m%d").date()  # filename date
        data_date_raw = datetime.strptime(line_items_step2[0][:2], "%d").date()    # date according to data
        data_time_raw = datetime.strptime(line_items_step2[0][2:], "%H%M").time() # time according to data

        if data_date_raw.day == filename_date_raw.day:
            today_utc_time = datetime.combine(filename_date_raw, data_time_raw)
        elif data_date_raw.day-1 == filename_date_raw.day:
            if data_time_raw.hour == 23:        # UTC 11pm the previous day
                print("Line of data in {} for station_id {} has data that is for the previous day".format(filename, station_id))
            else:       # If there are timestamps for the previous day but earlier than 11pm UTC
                print("Line of data in {} for station_id {} has observations outside of the 1 hour from the previous day".format(filename, station_id))
                continue
        else:        # NOTE: some are: HH-MM-SS ???? they definitely don't match up with what the file date is -- throwing out for now
            print("Line of data in {} for station_id {} has observations that do not match the date.".format(filename, station_id))

        utc_time = today_utc_time
        print(station_id, utc_time)

        sfcWind_dir_raw = line_items_step2[3] # wind direction, degrees (from true north)
        sfcWind_raw = line_items_step2[4]  # wind speed, miles per hour
        sfcWind_gust_raw = line_items_step2[5] # wind gust, miles per hour
        tas_raw = line_items_step2[6][:3]  # air temperature, degF

        # Notes: APRS weather specification comments http://www.aprs.org/aprs11/spec-wx.txt
        var_opts = ["r", "P", "p", "h", "b", "L", "l"] # and "s" if necessary
        var_labs = []
        var_idx = []
        for item in var_opts:
            if item in line:
                item_label = item + "_raw"
                item_idx = line_items_step2[6][3:].index(item)
                var_labs.append(item_label)
                var_idx.append(item_idx)
            else:
                continue
        var_list = dict(zip(var_labs, var_idx))

        if 'P_raw' in var_list:  # rainfall since midnight, hundredths of an inch
            pr_raw = line_items_step2[6][var_list['P_raw']+4:var_list['P_raw']+7]
        else:
            # print("Variable: Precipitation since midnight not provided")
            continue

        if 'r_raw' in var_list:  # rainfall in last 1 hour, hundredths of an inch
            r_raw = line_items_step2[6][var_list['r_raw']+4:var_list['r_raw']+7]
        else:
            # print("Variable: Rainfall in the last 1 hour not provided")
            continue

        if 'p_raw' in var_list:  # rainfall in last 24 hours, hundredths of an inch
            p_raw = line_items_step2[6][var_list['p_raw']+4:var_list['p_raw']+7]
        else:
            # print("Variable: Rainfall in last 24 hours not provided")
            continue

        if 'h_raw' in var_list:     # relative humidity, %
            hurs_raw = line_items_step2[6][var_list['h_raw']+4:var_list['h_raw']+6]
        else:
            # print("Variable: Relative Humidity not provided")
            continue

        if 'b_raw' in var_list:     # barometric pressure, convert to standardized air pressure, UNIT DEPENDENT DECIMAL PLACE, currently hundredths of hPa (mbar)
            ps_raw = line_items_step2[6][var_list['b_raw']+4:var_list['b_raw']+9]
        else:
            # print("Variable: Air pressure not provided")
            continue

        if 'L_raw' in var_list:     # solar radiation, w/m2 -- L is for values below 999, l is above 1000
            rsds_raw = line_items_step2[6][var_list['L_raw']+4:var_list['L_raw']+7]
        elif 'l_raw' in var_list:
            rsds_raw = line_items_step2[6][var_list['l_raw']+4:var_list['l_raw']+7]
        else:
            # print("Variable: Solar radiation not provided")
            continue

        # # Snowfall is listed in the documentation but is seemingly not actually provided
        # if 's_raw' in var_list:     # snowfall in inches in last 24 hours
        #     s_raw = line_items_step2[6][var_list['s_raw']+4:var_list['s_raw']+7]
        # else:
        #     # print("Variable: snowfall not provided")
        #     continue

        print(station_id, hardware, utc_time_raw, lat_clean, lon_clean, sfcWind_dir_raw, sfcWind_raw, sfcWind_gust_raw, tas_raw, hurs_raw, ps_clean, rsds_raw)
        print(count)
        if not line: break
    print("{} cleaned".format(filepath))

parse_cwop_sr(file)



## Step 2: Drop unnecessary variables
drop_vars = [sfcWind_gust_raw, pr_raw, P_raw]  # only drops vars that are not priority variables
drop_vars_exclusive = []    # drops all priority vars except solar radiation -- this is a merge with CWOP full data question


# # Write errors to a csv file
# filepath = "errors_cwop_sr_vars_{}.csv".format(end_api) # Sets path for error file
# # print(errors)
# with open(filepath, "w") as outfile:
#     write = csv.writer(outfile)
#     writer.writerow(errors.keys)
#     writer.writerows(zip(*errors.values()))
#
# os.chdir(homedir)
#
# return dropvars

## Step 3: Converts station metadata to standard format, with unique identifier
network_name = "CWOP"
network_id = network_name + "_" + station_id


## Step 5: Converts missing data to standard format (this needs to come first because of CWOP_SR missing data flags)
all_vars = [sfcWind_dir_raw, sfcWind_raw, sfcWind_gust_raw, tas_raw, pr_raw, P_raw, r_raw, p_raw, b_raw, rsds_raw, hurs_raw]
for item in all_vars:
    if item == "..." or "....":     # air pressure is 4 chars
        item = -998
    else:
        item = float(item)      # converts string text to float


## Step 4: Converts data metadata to standard format, and converts units into standard units if not provided in standard units
## Hardcoding for now
tas_clean = ((tas_raw - 32.) * (5/9))   # deg F to deg C
ps_clean = ((float(ps_raw) / 10))   # decimal hPa to hPa
p_clean = (p_raw * 0.254)   # hundredths of inch to mm



## Step 6: Tracks existing qa/qc flag for review
flag1 = 998 # long sequences of 998 is a solar radiation issue
flag2 = -998 # missing data -- intermediary flag, original flag is "..."



## Step 7: Merge files by station, and outputs cleaned variables as a single .nc file for an individual network
