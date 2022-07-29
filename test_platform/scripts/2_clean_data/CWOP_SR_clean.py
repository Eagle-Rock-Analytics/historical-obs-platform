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

# set envr variables
datadir = "/historical-obs/historical-obs-platform/test_platform/CWOP_SR/"

# Step 1: Read through variables, and calculates variables if not observed
# Note: Solar data starts at 11pm UTC in each file - each file comprises mix of data from the named day + last hour of previous day

filename = "L20200101.txt"
file = os.path.join(datadir, filename)
date = open(file)
# print(date.readline())  # Getting line of text to see format

def parse_cwop_sr(filepath):
    """
    Parses solar radiation data for CWOP.
    Paramters: filepath (str): filepath for file to be parsed
    Returns: data (xr.array): parsed data
    """
    data = []

    with open(filepath, 'r') as f:
        lines.append(f.readlines())

    # line = "11655> 011315z 3902.33N /08711.98W _260 /000 g000 t030 P000 h91 b10134 L007 ws31"
    # line = "CW3702> 011845z 3907.05N /10441.50W _000 /000 g000 t044 r000 p000 P000 h45 b10009 L158 .DsWLL"

    # Splits line based off of the following string options, case matters
    for line in lines:
        # remove any item from the "hardware" list to ensure it doesn't interfere with variable extraction
        hardware_opts = ["DsIP", "AmbientCWOP.com", "eMB51", ".WD 31", "ws31", ".DsWLL", "DsVP",
                            "eCumulusDsVP", "WeatherCatV312B34H31", "eMB50", ".weewx-4.5.1-Vantage"]
        for item in hardware_opts:
            if item in line::
                line = line[:-len(item)]

        line_items = re.split(r"[>z/_gt]", line) ## come back to this later
        station_id = line_items[0]
        utc_time = line_items[1] # convert to datetime object, decimal UTC hours

        # does this step need to first?
        if line_items[2][7] == "S":
            continue
        elif lat = line_items[2]
            #convert lat

        if line_items[3][8] == "E":
            continue
        elif lon = line_items[3]
            # convert lon

        sfcWind_dir = line_items[4] # wind direction, degrees (from true north)
        sfcWind = line_items[5] # wind speed, miles per hour OR KNOTS??
        sfcWind_gust = line_items[6] # wind gust, miles per hour
        tas = line_items[7][:3]  # air temperature, degF

        # this is where additional variables tend to get added in, after air_temp
        # in better code, this would be automated as:
        # if any_character in a string of characters equals desired_character:
        # next 2-4 characters are saved to a particular variable
        # elif ... same thing but with a new character to look for

        if line_items[7][3:] == "r":    # idk
            r_idx = line_items[7][3:].index("r")
            rvar = line_items[7][r_idx+4:r_idx+7]
        elif line_items[7][3:] == "P":  # precipitation, hundredths of an inch
            P_idx = line_items[7][3:].index("P")
            pr = line_items[7][P_idx+4:P_idx+7]
        elif line_items[7][3:] == "p":  # idk
            p_idx = line_items[7][3:].index("p")
            pvar = line_items[7][p_idx+4:p_idx+7]
        elif line_items[7][3:] == "h":  # relative humidity, %
            h_idx = line_items[7][3:].index("h")
            hurs = line_items[7][h_idx+4:h_idx+6]
        elif line_items[7][3:] == "b":  # barometric pressure, convert to standardized air pressure, UNIT DEPENDENT DECIMAL PLACE, currently hundredths of hPa
            b_idx = line_items[7][3:].index("b")
            ps = line_items[7][b_idx+4:b_idx+8]
        elif line_items[7][3:] == "L": # solar radiation, w/m2 -- L is for values below 999 ## known flag is stretches of 998
            L_idx = line_items[7][3:].index("L")
            rsds = line_items[7][L_idx+4:L_idx+7]
        elif line_items[7][3:] == "l": # solar radiation, w/m2 -- l is for values above 1000 ## known flag is stretches of 998
            l_idx = line_items[7][3:].index("l")
            rsds = line_items[7][l_idx+4:l_idx+7]


parse_cwop_sr(date)



# Step 2: Drop unnecessary variables
drop_vars = []



# Step 3: Converts station metadata to standard format, with unique identifier
network_name = "CWOP"
station_id = network_name + "_" + network_id

# Step 4: Converts data metadata to standard format, and converts units into standard units if not provided in standard units.

# unit conversions
tas = (5/9) * (tas - 32) + 273.15 # convert to K
tas.units = "K"



# Step 5: Converts missing data to standard format



# Step 6: Tracks existing qa/qc flag for review
flag1 = 998 # long sequences of 998 is a solar radiation issue
flag2 = '...' # missing data



# Step 7: Merge files by station, and outputs cleaned variables as a single .nc file for an individual network
