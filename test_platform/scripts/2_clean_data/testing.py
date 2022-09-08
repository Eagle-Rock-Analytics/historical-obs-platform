import os
from datetime import datetime, timezone, timedelta
import xarray as xr
import re
import pandas as pd
from pathlib import Path
import csv
import calc_clean
import string
import numpy as np


homedir = os.getcwd() # Get current working directory.
if "historical-obs-platform" in homedir: # If git folder in path
    homedir = homedir[0:homedir.index("historical-obs-platform")]+"historical-obs-platform" # Set path to top folder.
    os.chdir(homedir) # Change directory.
else:
    print("Error: Set current working directory to the git repository or a subfolder, and then rerun script.")
    exit()

raw_datadir = homedir + "/test_platform/data/1_raw_wx/CWOP_SR/"
clean_datadir = homedir + "/test_platform/data/2_clean_wx/CWOP_SR/"

f_ = os.listdir(raw_datadir)
f_raw = list(filter(lambda f: f.endswith(".txt"), f_))

filename = "L20200515.txt"
file = os.path.join(raw_datadir, filename) # pathway to file

os.chdir(clean_datadir)

try:
    t, m, bbox = calc_clean.get_wecc_poly(wecc_terr, wecc_mar)
    lonmin, lonmax = float(bbox['minx']), float(bbox['maxx'])
    latmin, latmax = float(bbox['miny']), float(bbox['maxy'])
except:
    lonmin, lonmax = -139.047795, -102.03721
    latmin, latmax = 30.142739, 60.003861

gd = 0
good_stns = []
all_stns = []

f = open(file, "r")
for line in f:
    line_items_step1 = re.split(r"[>zh/_]", line) # Strip off station_id first
    stn_id = line_items_step1[0]
    all_stns.append(stn_id)
    utc_raw = line_items_step1[1] # [6] is either z or h
    lat_raw = line_items_step1[2] # [15] is /
    lon_raw = line_items_step1[3] # [25] is _

    # Longitude is stable
    if lon_raw[-1] == "W" or lon_raw[-1] == "w": # This also catches some "bad coded" missing lat-lon coords
        if lon_raw[3] == "." and lon_raw[6] == ".": # XX.XX.XX format -- separate DMS format
            lon_clean = calc_clean._lon_dms_to_dd(lon_raw[:-1])
        elif lon_raw[1] == "." and lon_raw[17] == "+": # X.Xe+00X notation.... apparently just one station - but it is in WECC
            lon_raw = lon_raw[:-1]
            _deg = float(lon_raw[:3]) * 10
            _min = float(lon_raw[22:25]) + float(lon_raw[-3:])
            lon_clean = -1 * (_deg + (_min)/60)
        elif lon_raw[5] != ".": # XX.XXXX format -- already in Dd but some obs are badly formatted (XX instead of XXX for lon)
            lon_clean = -1 * float(lon_raw[:-1])
            # print(stn_id, lon_raw, lon_clean)
        else: # LORAN format (DDMM.mm) -- most obs will fall in this category
            lon_clean = calc_clean._lon_DMm_to_Dd(lon_raw[:-1])
    else:
        continue # Skips if it is the wrong hemisphere

    if lon_clean < lonmax and lon_clean > lonmin:
        lon_clean = lon_clean # inside WECC
    elif lon_clean > lonmax:
        continue # east of WECC
    else:
        continue # west of WECC


    # Latitude is stable
    if len(lat_raw) > 4:    # Specifically catches some "bad coded" missing lat-lon coords
        if lat_raw[0] != "-" and (lat_raw[-1] == "N" or lat_raw[-1] == "n"): # More catching of bad coded lat-lon coords
            if lat_raw[2] == "." and lat_raw[5] == ".": # XX.XX.XX format -- separate DMS format
                lat_clean = calc_clean._lat_dms_to_dd(lat_raw[:-1])
            elif lat_raw[1] == "." and lat_raw[17] == "+":  # X.Xe+00X notation.... apparently just one station - but it is in WECC
                lat_raw = lat_raw[:-1]
                _deg = float(lat_raw[:3]) * 10
                _min = float(lat_raw[23:26]) + float(lat_raw[-3:])
                lat_clean = _deg + _min/60
            elif lat_raw[2] == "." and lat_raw[5] != ".": # XX.XXXX format -- already in Dd but some obs are badly formatted
                lat_clean = float(lat_raw[:-1])
            elif lat_raw[2] == ",": # europeans
                lat_raw = lat_raw[:-1]
                lat_clean = float(lat_raw[:2]) + float(lat_raw[3:])/100
            else:   # LORAN format (DDMM.mm) -- most obs will fall in this category
                lat_clean = calc_clean._lat_DMm_to_Dd(lat_raw[:-1])
        else:
            continue # wrong hemisphere
    else:
        continue

    if lat_clean < latmax and lat_clean > latmin:
        lat_clean = lat_clean
        gd += 1
        good_stns.append(stn_id)
    else:
        continue # south or north of WECC

    good_data = gd

#-------------------------------------------------------------------------------
print('Usable data: ', good_data)

fp = open(file, "r")
for count, line in enumerate(fp):
    pass
print('Total Lines of data: ', count + 1)
print('Percentage usable: ', good_data/(count+1))

g_res = np.array(good_stns)
good_cwop_stns_n = len(np.unique(g_res))

res = np.array(all_stns)
all_cwop_stns_n = len(np.unique(res))

print("Usable # of stations: ", good_cwop_stns_n)
print("Total # of stations : ", all_cwop_stns_n)
print("Percentage of usable stations: ", good_cwop_stns_n / all_cwop_stns_n)

#-------------------------------------------------------------------------------

os.chdir(homedir + "/test_platform/scripts/2_clean_data/")
cwop_data = pd.read_csv("cwop_stations.csv", usecols=["STID"])

cwop_data_arr = np.array(cwop_data)

overlap = np.isin(cwop_data_arr, np.unique(g_res), assume_unique=True)
stn_ct = 0
for i in overlap:
    if i == True:
        # print('there is a station in CWOP_SR in CWOP')
        stn_ct += 1
    else:
        # print('there is no station overlap between CWOP_SR and CWOP')
        continue
# print(overlap)
print('Number of CWOP_SR stations in CWOP: ', stn_ct)
print('Those stations are: ', cwop_data_arr[overlap])
