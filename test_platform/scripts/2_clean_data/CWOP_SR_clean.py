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

# testing using 2020 and 2021

# Step 0: Environment set-up
# Import libraries
import os
from datetime import datetime, timezone
import xarray as xr
import re
import pandas as pd
# from .thermo_helpers import f_to_c

# set envr variables
datadir = "/Users/victoriaford/Desktop/historical-obs/historical-obs-platform/test_platform/CWOP_SR/"

# Step 1: Read through variables, and calculates variables if not observed
# Note: Solar data starts at 11pm UTC in each file - each file comprises mix of data from the named day + last hour of previous day

# Hardcoding boundaries in for now (can use WECC poly in future)
lonmin, lonmax = -139.047795, -102.03721
latmin, latmax = 30.142739, 60.003861

filename = "L20200101.txt"
file = os.path.join(datadir, filename) # pathway to file
date = open(file)
# print(date.readline())  # Getting line of text to see format

def parse_cwop_sr(filepath):
    """
    Parses solar radiation data for CWOP.
    Paramters: filepath (str): filepath for file to be parsed
    Returns: data (xr.array): parsed data
    """
    lines = []

    date = open(filepath)

    f = open(filepath, 'r')
    count = 0

    while True:
        count += 1
        line = f.readline()

        line_items = re.split(r"[>z/_gt]", line)
        station_id = line_items[0] ## FLAGGING: CW8876 has duplicate station_id
        utc_time = line_items[1]

        # Lat-lon conversion: stip hemisphere designator, convert to decimal degrees, removes data outside WECC region of interest
        if line_items[2][7] == "N":
            lat_raw = float(line_items[2][:-1]) # Strips the N char off and converts to float
            deg = float(line_items[2][:2])
            min_ = float(line_items[2][2:4])
            sec_ = float(line_items[2][5:7])
            lat_clean = deg + min_/60 + sec_/3600
            if lat_clean < latmin: # getting rid of latitude locations outside WECC region (hard coded)
                continue
            elif lat_clean > latmax:
                continue
            else:
                lat_clean = lat_clean
        elif line_items[2][7] == "S":
            continue
        else:   # primarily error handling in case a line does not have N or S provided -- is this the case or can we get rid of?
            continue
            print("Error: this latitude does not exist/incorrect or missing hemisphere designator.")

        if line_items[3][8] == "W":
            lon_raw = float(line_items[3][:-1]) # Strips the W char off and converts to float
            deg = float(line_items[3][:3])
            min_ = float(line_items[3][3:5])
            sec_ = float(line_items[3][6:8])
            lon_clean = (deg + min_/60 + sec_/3600) * -1 # converting to negative values for standard western hemi reporting
            if lon_clean < lonmin:
                continue
            elif lon_clean > lonmax:
                continue
            else:
                lon_clean = lon_clean
        elif line_items[3][8] == "E":
            continue
        else:   # primarily error handling in case a line does not have E/W provided -- is this the case or can we get rid of?
            continue
            print("Error: this longitude does not exist/incorrect or missing hemisphere designator.")

        ## NOTE: CWOP_SR uses "..." as a missing data flag for QA/QC
        if line_items[4] == "...":
            sfcWind_dir = "..." # Keeping as string for now
        else:
            sfcWind_dir = float(line_items[4]) # wind direction, degrees (from true north

        ## Check wind speed unit consistency -- are stations using mph or knots?
        if line_items[5] == "...":
            sfcWind = "..."
        else:
            sfcWind = float(line_items[5])  # wind speed, miles per hour OR KNOTS??

        if line_items[6] == "...":
            sfcWind_gust = "..."
        else:
            sfcWind_gust = float(line_items[6]) # wind gust, miles per hour

        if line_items[7][:3] == "...":
            tas = "..."
        else:
            tas_raw = float(line_items[7][:3])  # air temperature, degF

        # print(station_id, utc_time, lat_clean, lon_clean, sfcWind_dir, sfcWind, sfcWind_gust, tas_clean) # Testing check

        # r is ....
        # p is .... (rainfall per day?)
        var_opts = ["r", "P", "p", "h", "b", "L", "l"]
        for item in var_opts:
            if item in line:
                item_idx = line_items[7][3:].index(item)
                item_label = item + "_raw"
                if P_raw == True:
                    pr_raw = line_items[7][P_idx+4:P_idx+7] # precipitation, hundredths of an inch
                elif r_raw == True:
                    r_raw = line_items[7][r_idx+4:r_idx+7]  # idk
                elif p_raw == True:
                    p_raw = line_items[7][p_idx+4:p_idx+7]  # idk
                elif h_raw == True:
                    hurs_raw = line_items[7][h_idx+4:h_idx+6]  # relative humidity, %
                elif b_raw == True:
                    ps_raw = line_items[7][b_idx+4:b_idx+8]  # barometric pressure, convert to standardized air pressure, UNIT DEPENDENT DECIMAL PLACE, currently hundredths of hPa
                elif L_raw == True:
                    rsds_raw = line_items[7][L_idx+4:L_idx+7]  # solar radiation, w/m2 -- L is for values below 999
                elif l_raw == True:
                    rsds_raw = line_items[7][l_idx+4:l_idx+7]  # solar radiation, w/m2 -- L is for values below 999
                else:
                    print("Error -- unknown variable flag present in data, please confirm")


        # hardware_opts = ["DsIP", "AmbientCWOP.com", "eMB51", ".WD 31", "ws31", ".DsWLL", "DsVP", "DsIP", "WD", "eCumulus",
        #                 "eCumulusDsVP", "WeatherCatV312B34H31", "eMB50", ".weewx-4.5.1-Vantage", "eMB", "weewx", "eMH", "WeatherCat"
        #                 "eCumulusD", "wview", "Vantage", "Davis", "VP"] ## technical report indicates there are more than this -- need to flag/automate
        #
        # for item in hardware_opts:
        #     if item in line:
        #         line_r = line[:-len(item)]
        #         print(line_r)


        if not line:
            break
        # print("Line: {}".format(count, line.strip()))
    # lines.append(f.readlines())

    # line = "11655>011315z3902.33N/08711.98W_260/000g000t030P000h91b10134L007ws31"
    # line = "CW3702>011845z3907.05N/10441.50W_000/000g000t044r000p000P000h45b10009L158.DsWLL"

    # Splits line based off of the following string options, case matters
    # for line in lines:
    #     # remove any item from the "hardware" list to ensure it doesn't interfere with variable extraction
    #     hardware_opts = ["DsIP", "AmbientCWOP.com", "eMB51", ".WD 31", "ws31", ".DsWLL", "DsVP",
    #                         "eCumulusDsVP", "WeatherCatV312B34H31", "eMB50", ".weewx-4.5.1-Vantage"] ## technical report indicates there are more than this -- need to flag/automate
    #     for item in hardware_opts:
    #         if item in line:
    #             line = line[:-len(item)]
    #             print(line)
    #         else:
    #             continue


parse_cwop_sr(file)

######
#
# # Step 2: Drop unnecessary variables
drop_vars = [sfcWind_gust]  # only drops vars that are not priority variables

drop_vars_exclusive = []    # drops all priority vars except solar radiation -- this is a merge with CWOP full data question
#
#
#
# Step 3: Converts station metadata to standard format, with unique identifier
network_name = "CWOP"
station_id = network_name + "_" + network_id

# Step 4: Converts data metadata to standard format, and converts units into standard units if not provided in standard units.
# unit conversions
# tas_clean = f_to_c(tas_raw)
tas_clean = ((tas_raw - 32.) * (5/9))   # hard coding for now
# print(tas_raw, tas_clean)   # conversion check


# Step 5: Converts missing data to standard format



# Step 6: Tracks existing qa/qc flag for review
flag1 = 998 # long sequences of 998 is a solar radiation issue
flag2 = '...' # missing data



# Step 7: Merge files by station, and outputs cleaned variables as a single .nc file for an individual network
