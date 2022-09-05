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

## Testing using 2020 and 2021

## Step 0: Environment set-up
## Import libraries
import os
from datetime import datetime, timezone, timedelta
import xarray as xr
import re
import pandas as pd
from pathlib import Path
import csv
import calc_clean

## set envr variables
homedir = os.getcwd() # Get current working directory.
if "historical-obs-platform" in homedir: # If git folder in path
    homedir = homedir[0:homedir.index("historical-obs-platform")]+"historical-obs-platform" # Set path to top folder.
    os.chdir(homedir) # Change directory.
else:
    print("Error: Set current working directory to the git repository or a subfolder, and then rerun script.")
    exit()

raw_datadir = homedir + "/test_platform/data/1_raw_wx/CWOP_SR/"
clean_datadir = homedir + "/test_platform/data/2_clean_wx/CWOP_SR/"
wecc_terr = homedir + "/test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = homedir + "/test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"


## Set up directory to save files, if it doesn't already exist.
try:
    os.mkdir(clean_datadir) # Make folder to save cleaned data
    print("Directory for {} created".format(re.split("/", raw_datadir)[-2]))
except:
    print("Directory for {} exists".format(re.split("/", raw_datadir)[-2]))
    pass    # Pass if folder already exists


## Step 1: Read through variables, and calculates variables if not observed
filename = "L20200102.txt"
file = os.path.join(raw_datadir, filename) # pathway to file

#------------------------------------------------------------------------------------------------------------
## Get CWOP SOLAR RAD STATION ELEVATIONS
def get_cwop_elevs(filepath):
    # calculate elevation -- from CWOP_SR they use a website to determine elevations -- url api?
    print("Testing: Grabbing elevations for stations")

get_cwop_elevs(file)

#------------------------------------------------------------------------------------------------------------
## Parse through variables
def parse_cwop_sr(filepath):
    """
    Parses solar radiation data for CWOP.
    Paramters: filepath (str): filepath for file to be parsed
    Returns: data (xr.array): parsed raw data
    """

    os.chdir(clean_datadir) # Switch to cleaning directory
    try:
        t, m, bbox = calc_clean.get_wecc_poly(wecc_terr, wecc_mar)
        lonmin, lonmax = float(bbox['minx']), float(bbox['maxx'])
        latmin, latmax = float(bbox['miny']), float(bbox['maxy'])
    except:
        lonmin, lonmax = -139.047795, -102.03721
        latmin, latmax = 30.142739, 60.003861

    ## Set up csv to record any files not cleaned and reason why
    errors = {'File':[], 'Time':[], 'Error':[]}
    end_api = datetime.now().strftime("%Y%m%d%H%M") # sets end time to be current time at beginning of download

    os.chdir(raw_datadir) # Change directory to where raw data files are stored
    files = os.listdir()
    files = list(filter(lambda f: f.endswith(".txt"), files))
    # print(files)

    ## Obtaining station elevations  -- IN PROGRESS
    try:
        # url = "https://www.ndbc.noaa.gov/bmanht.shtml"
        elevs = get_cwop_elevs()
    except Exception as e:
        print('Testing: Reading elevations for stations')
        # continue

    for i in ['L20200102']: #, 'L20200629', 'L20201231']:  # TESTING
        print("Cleaning: {}".format(i))
        ds_stat = None

        f = open(filepath, "r")
        # line = f.readline()

        try:
            for line in f:
                # Strips off the leading station_id, not of equal character length
                line_items_step1 = re.split(r">", line) # Strip off station_id first
                station_id = line_items_step1[0]
                ## Removing data outside WECC region first before any other conversions
                ## Lat-lon conversion: strip hemisphere designator, convert to decimal degrees, removes data outside WECC region of interest
                ## [0-6] is utc time -- can't strip first because of bad-data flags in lat-lon
                lat_raw = line_items_step1[1][7:15] # [15] is /
                lon_raw = line_items_step1[1][16:25] # [25] is _
                try:
                    if lat_raw[-1:] == "N" or lat_raw[-1:] == "n":
                        if calc_clean.uses_other_chars(lat_raw) == False: # data is good
                            if lat_raw[2] != ".":
                                lat_clean = calc_clean._lat_DMm_to_Dd(lat_raw[:-1])
                            else:
                                if lat_raw[5] == ".":
                                    lat_clean = lat_raw[:-1]
                                else:
                                    lat_clean = calc_clean._lat_dms_to_dd(lat_raw[:-1])
                        if lat_clean < latmin or lat_clean > latmax: # getting rid of latitude locations outside WECC region
                            errors['File'].append(filename)
                            errors['Time'].append(end_api)
                            errors['Error'].append("Line of data not in WECC. Lat: {}".format(lat_clean))
                            continue
                        else:
                            lat_clean = lat_clean
                    else:
                        errors['File'].append(filename)
                        errors['Time'].append(end_api)
                        errors['Error'].append("Line of data not in WECC. Lat: {}".format(lat_clean))
                        continue

                    if lon_raw[-1:] == "W" or lon_raw[-1:] == "w":
                        if calc_clean.uses_other_chars(lon_raw) == False: # data is good
                            if lon_raw[3] == ".":
                                lon_clean = calc_clean._lon_dms_to_dd(lon_raw[:-1])
                            else:
                                lon_clean = calc_clean._lon_DMm_to_Dd(lon_raw[:-1])
                        else:
                            print('Bad Data Coding for Longitude') # bad data
                        # if lon_clean < lonmin or lon_clean > lonmax:  # get_wecc_poly
                        # if lon_clean > -139.047795 and lon_clean < -102.03721: # hard coding, since get_wecc_poly not behaving
                        #     lon_clean = lon_clean
                        # else:
                        #     print('not in wecc poly')
                        #     errors['File'].append(filename)
                        #     errors['Time'].append(end_api)
                        #     errors['Error'].append("Line of data not in WECC. Lon: {}".format(lon_clean))
                        #     continue
                    else:
                        errors['File'].append(filename)
                        errors['Time'].append(end_api)
                        errors['Error'].append("Line of data not in WECC. Lon: {}".format(lon_clean))
                        continue
                except Exception as e:
                    print('something is breaking here')
                    continue # Continue to next branch

                # Writing to file after geographic subsetting
                stn_filepath = clean_datadir + "cwop_solarrad_stations.csv"
                with open(stn_filepath, "w") as outfile:
                    writer = csv.writer(outfile)
                    writer.writerow([station_id]) # Would be great if was station_id, first date, last date of coverage

                ## Note: Solar data starts at 11pm UTC in each file - each file comprises mix of data from the named day + last hour of previous day
                ## Keep the data from previous day, because it will get merged by station anyways
                filename_date_raw = datetime.strptime(filename[1:-4].strip(), "%Y%m%d").date()  # filename date
                yesterday_date_raw = filename_date_raw + timedelta(days=-1) # previous day date - for 11pm UTC tracking
                utc_time_idx = line_items_step1[1][:6]
                data_date_raw = datetime.strptime(utc_time_idx[:2], "%d").date()    # date according to data
                data_time_raw = datetime.strptime(utc_time_idx[2:], "%H%M").time() # time according to data

                try:
                    if data_date_raw.day == filename_date_raw.day:
                        utc_time_clean = datetime.combine(filename_date_raw, data_time_raw)
                    elif data_date_raw.day == yesterday_date_raw.day:
                        utc_time_clean = datetime.combine(yesterday_date_raw, data_time_raw)
                except:        # Bad data date (day does not match filename or the previous day)
                    errors['File'].append(filename)
                    errors['Time'].append(end_api)
                    errors['Error'].append("Line of data in {} for station_id {} has observations that do not match the date.".format(filename, station_id))
                    continue

                ## Strip Hardware off here -- because of unequal variable counts, and duplicating letters
                rad_char = "L"
                rad_pos = line_items_step1[1][25:].lower().find(rad_char.lower())   # case insensitive
                hardware = line_items_step1[1][25+rad_pos+4:]

                line_items_step2 = line_items_step1[1][25:25+rad_pos+4] # Strips off Hardware
                line_items_step3 = re.split(r"[/gt]", line_items_step2) # starts at the _ before wind direction
                sfcWind_dir_raw = line_items_step3[0][1:] # wind direction, degrees (from true north)
                sfcWind_raw = line_items_step3[1]  # wind speed, miles per hour
                sfcWind_gust_raw = line_items_step3[2] # wind gust, miles per hour
                tas_raw = line_items_step3[3][:3]  # air temperature, degF

                # Note: after 't', order of variables can vary
                # Note: APRS weather specification comments http://www.aprs.org/aprs11/spec-wx.txt
                var_opts = ["r", "P", "p", "h", "b", "L", "l"]  # and "s" if necessary
                var_labs = []
                var_idx = []
                for item in var_opts:
                    if item in line:
                        item_label = item + "_raw"
                        item_idx = line_items_step3[-1][3:].index(item)
                        var_labs.append(item_label)
                        var_idx.append(item_idx)
                    else:
                        continue
                var_list = dict(zip(var_labs, var_idx))
                print(var_list)     ## testing

                try:
                    if 'P_raw' in var_list:  # rainfall since midnight, hundredths of an inch -- drop
                        pr_raw = line_items_step3[-1][var_list['P_raw']+4:var_list['P_raw']+7]
                    else:
                        errors['File'].append(filename)
                        errors['Time'].append(end_api)
                        errors['Error'].append("Line of data in {} for station_id {} does not have observations for pr_raw.".format(filename, station_id))
                        continue

                    if 'r_raw' in var_list:  # rainfall in last 1 hour, hundredths of an inch -- keep
                        r_raw = line_items_step3[-1][var_list['r_raw']+4:var_list['r_raw']+7]
                    else:
                        errors['File'].append(filename)
                        errors['Time'].append(end_api)
                        errors['Error'].append("Line of data in {} for station_id {} does not have observations for r_raw.".format(filename, station_id))
                        continue

                    if 'p_raw' in var_list:  # rainfall in last 24 hours, hundredths of an inch -- drop
                        p_raw = line_items_step3[-1][var_list['p_raw']+4:var_list['p_raw']+7]
                    else:
                        errors['File'].append(filename)
                        errors['Time'].append(end_api)
                        errors['Error'].append("Line of data in {} for station_id {} does not have observations for p_raw.".format(filename, station_id))
                        continue

                    if 'h_raw' in var_list:     # relative humidity, %
                        hurs_raw = line_items_step3[-1][var_list['h_raw']+4:var_list['h_raw']+6]
                    else:
                        errors['File'].append(filename)
                        errors['Time'].append(end_api)
                        errors['Error'].append("Line of data in {} for station_id {} does not have observations for hurs_raw.".format(filename, station_id))
                        continue

                    if 'b_raw' in var_list:     # barometric pressure, convert to standardized air pressure, UNIT DEPENDENT DECIMAL PLACE, currently hundredths of hPa (mbar)
                        ps_raw = line_items_step3[-1][var_list['b_raw']+4:var_list['b_raw']+9]
                    else:
                        errors['File'].append(filename)
                        errors['Time'].append(end_api)
                        errors['Error'].append("Line of data in {} for station_id {} does not have observations for ps_raw.".format(filename, station_id))
                        continue

                    if 'L_raw' in var_list:     # solar radiation, w/m2 -- L is for values below 999, l is above 1000
                        L_rsds_raw = line_items_step3[-1][var_list['L_raw']+4:var_list['L_raw']+7]
                    else:
                        # print('No L')
                        errors['File'].append(filename)
                        errors['Time'].append(end_api)
                        errors['Error'].append("Line of data in {} for station_id {} does not have observations for rsds_raw.".format(filename, station_id))
                        continue

                    if 'l_raw' in var_list:
                        l_rsds_raw = line_items_step3[-1][var_list['l_raw']+4:var_list['l_raw']+7]
                    else:
                        # print('No l')
                        errors['File'].append(filename)
                        errors['Time'].append(end_api)
                        errors['Error'].append("Line of data in {} for station_id {} does not have observations for rsds_raw.".format(filename, station_id))
                        continue

                except Exception as e:
                    errors['File'].append(filename)
                    errors['Time'].append(end_api)
                    errors['Error'].append("Line of data in {} for station_id {} does not have observations for {}.".format(filename, station_id, utc_time_clean))
                    continue

                ## Snowfall is listed in the documentation but is seemingly not actually provided
                # if 's_raw' in var_list:     # snowfall in inches in last 24 hours
                #     s_raw = line_items_step3[-1][var_list['s_raw']+4:var_list['s_raw']+7]


            ## Write files to a new netcdf
            # desired_order = ['ps', 'tas', 'tdps', 'pr', 'hurs', 'rsds', 'sfcWind', 'sfcWind_dir']
            # try:
            #     filename = id + ".nc"
            #     filepath = clean_datadir + filename
            #     ds_stat.to_netcdf(path = filepath)
            #     print("Saving {} with dims {}".format(filename, ds_stat.dims))
            #     ds.close()
            # except Exception as e:
            #     print(filename, e)
            #     errors['File'].append(filename)
            #     errors['Time'].append(end_api)
            #     errors['Error'].append(e)
            #     continue

            f.close()
            print("{} cleaned".format(filename))
            # print("{} lines of data within WECC poly region cleaned".format()) ## Great to provide for sense of how much data is in each text file

        except Exception as e:
            print(file, e) # testing only purposes
            errors['File'].append(file)
            errors['Time'].append(end_api)
            errors['Error'].append(e)
            continue

    # Write errors to a csv for tracking
    filepath = clean_datadir + "errors_cwop_solarrad_{}.csv".format(end_api)
    with open(filepath, "w") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(errors.keys())
        writer.writerows(zip(*errors.values()))

parse_cwop_sr(file)


        # #----------------------------------
        #
        #     id = "CWOP_SR_"+ds.attrs["id"]
        #     ds['original_id'] = ds.attrs['id'] # Keep original ID as variable
        #
        #     ## Set attributes for dimensional data.
        #     ## Time
        #     ds['time'].attrs['long_name'] = "time"
        #     ds['time'].attrs['standard_name'] = "time"
        #     ds['time'].attrs['comment'] = "In UTC."
        #
        #     ## Station ID
        #     ds['station'].attrs['long_name'] = "station_id"
        #     ds['station'].attrs['comment'] = "Unique ID created by Eagle Rock Analytics. Includes network name appended to original unique station ID provided by network."
        #
        #     ## Latitude
        #     ds['lat'].attrs['long_name'] = "latitude"
        #     ds['lat'].attrs['standard_name'] = "latitude"
        #     ds['lat'].attrs['units'] = "degrees_north"
        #
        #     ## Longitude
        #     ds['lon'].attrs['long_name'] = "longitude"
        #     ds['lon'].attrs['standard_name'] = "longitude"
        #     ds['lon'].attrs['units'] = "degrees_west" ### CHECK
        #
        #
        #     ## Step 4: Convert dataset metadata in standard format -- CF compliance - to be finalized, overwrite existing metadata.
        #     ds = ds.assign_attrs(title = "CWOP Solar Radiation cleaned")
        #     ds = ds.assign_attrs(institution = "Eagle Rock Analytics / Cal Adapt")
        #     ds = ds.assign_attrs(source = "")
        #     ds = ds.assign_attrs(history = "CWOP_solarrad_clean.py script run on {} UTC".format(timestamp))
        #     ds = ds.assign_attrs(comment = "Intermediate data product: may not have been subject to any cleaning or QA/QC processing")
        #     ds = ds.assign_attrs(license = "")
        #     ds = ds.assign_attrs(citation = "")
        #     ds = ds.assign_attrs(disclaimer = "This document was prepared as a result of work sponsored by the California Energy Commission (PIR-19-006). It does not necessarily represent the views of the Energy Commission, its employees, or the State of California. Neither the Commission, the State of California, nor the Commission's employees, contractors, or subcontractors makes any warranty, express or implied, or assumes any legal liability for the information in this document; nor does any party represent that the use of this information will not infringe upon privately owned rights. This document has not been approved or disapproved by the Commission, nor has the Commission passed upon the accuracy of the information in this document.")
        #
        #
        #     ## Step 4: Convert data metadata to standard format, and converts units into standard units if not provided in standard units
        #     ## If not observed, calculate derived primary variables if possible.
        #
        #     # tas: air surface temperature (fahrenheit)
        #     if "air_temperature" in ds.keys():
        #         ds = ds.rename({'air_temperature': 'tas'})
        #         ds['tas_raw'] = ds['tas']
        #         try:
        #             ds['tas'] = calc_clean._unit_degF_to_K(ds['tas'])
        #         except:
        #             print('tas: calc_clean.py not working.')
        #         ds['tas'] = (5/9) * (ds['tas'] - 32.) + 273.15  # backup
        #         ds['tas'].attrs['ancillary_variables'] = "tas_raw"
        #         ds['tas'].attrs['comment'] = "Converted from Fahrenheit."
        #         ds['tas'].attrs['long_name'] = "air_temperature"
        #         ds['tas'].attrs['standard_name'] = "air_temperature"
        #         ds['tas'].attrs['units'] = "degree_Kelvin"
        #
        #     # ps: surface air pressure (decimal hectopascals)
        #     if "air_pressure" in ds.keys():
        #         ds = ds.rename({'air_pressure': 'ps'})
        #         ds['ps_raw'] = ds['ps']
        #         try:
        #             ds['ps'] = calc_clean._unit_pres_hpa_to_pa(ds['ps']/10)
        #         except:
        #             print('ps: calc_clean.py not working')
        #         ds['ps'].attrs['ancillary_variables'] = "ps_raw"
        #         ds['ps'].attrs['comment'] = "Converted from decimal hectopascals."
        #         ds['ps'].attrs['long_name'] = "station_air_pressure"
        #         ds['ps'].attrs['standard_name'] = "air_pressure"
        #         ds['ps'].attrs['units'] = "Pa"
        #
        #     # pr: precipitation (rainfall in last 1 hour, hundredths of an inch)
        #     if "precipitation" in ds.keys():
        #         ds = ds.rename({'precipitation': "pr"})
        #         ds['pr_raw'] = ds['pr']
        #         try:
        #             ds['pr'] = calc_clean._unit_precip_in_to_mm(ds['pr'] * 100.)
        #         except:
        #             print('pr: calc_clean.py not working')
        #         ds['pr'].attrs['ancillary_variables'] = "pr_raw"
        #         ds['pr'].attrs['comment'] = "Converted from hundredths of an inch."
        #         ds['pr'].attrs['long_name'] = "precipitation"
        #         ds['pr'].attrs['standard_name'] = "precipitation"
        #         ds['pr'].attrs['units'] = "mm"
        #
        #     # relative humiidity (%)
        #     if "relative_humidity" in ds.keys():
        #         ds = ds.rename({'relative_humidity': 'hurs'})
        #         ds['hurs_raw'] = ds['hurs']
        #         ds['hurs'].attrs['long_name'] = "relative_humdidity"
        #         ds['hurs'].attrs['standard_name'] = "relative_humidity"
        #         ds['hurs'].attrs['units'] = "percent"
        #
        #     # dew point temperature (kelvin)
        #     if "dewpoint_temperature" in ds.keys():
        #         ds = ds.rename({'dewpoint_temperature': 'tdps'})
        #         ds['tdps_raw'] = ds['tdps']
        #         try:
        #             ds['tdps'] = calc_clean._calc_dewpointtemp_opt1(ds['tas'], ds['hurs'])
        #         except:
        #             print('tdps: calc_clean.py not working.')
        #         ds['tdps'].attrs['ancillary_variables'] = "tdps_raw"
        #         ds['tdps'].attrs['comment'] = "Calculated using primary method with air temperature and relative humidity"
        #         ds['tdps'].attrs['long_name'] = "dewpoint_temperature"
        #         ds['tdps'].attrs['standard_name'] = "dewpoint_temperature"
        #         ds['tdps'].attrs['units'] = "degree_Kelvin"

        #     if not line:
        #         break
        #
        # # Writes errors to csv
        # error_filepath = clean_datadir + "errors_cwop_solarrad_{}.csv".format(end_api)
        # with open(error_filepath, "w") as outfile:
        #     writer = csv.writer(error_filepath)
        #     writer.writerow(errors.keys())
        #     writer.writerows(zip(*errors.values()))
# print("{} cleaned".format(filename))
# parse_cwop_sr(file)

#------------------------------------------------------------------------------------------------------------
# def get_ids(raw_datadir):
#     """
#     Pulls in station list, and returns list of station IDs.
#     TO DO: compare CWOP and CWOP SOLAR RAD stations for overlapping cwop_stations
#     """
#     print(station_id)
#
#     filepath = "cwop_solarrad_stations.csv"
#     with open(filepath, "w") as outfile:
#         writer = csv.writer(outfile)
#         writer.writerow(station_id)
#
#     ## os.chdir(clean_datadir)
#     ## cwop_sr = pd.read_csv("cwop_solarrad_stations.csv")
#     ## cwop = pd.read_csv("cwop_stations.csv")
#     ##
#     ## for item in cwop_sr:
#     ##     if item not found in cwop
#     ##     print('{} station not found in cwop station list'.format(item))
#
# test = get_ids(raw_datadir)
# print(test) ###

#------------------------------------------------------------------------------------------------------------

# ## Step 2: Drop unnecessary variables
# drop_vars = [sfcWind_gust_raw, pr_raw, P_raw]  # only drops vars that are not priority variables
# drop_vars_exclusive = []    # drops all priority vars except solar radiation -- this is a merge with CWOP full data question
#
#
# ## Step 3: Converts station metadata to standard format, with unique identifier
# network_name = "CWOP_SR"
# network_id = network_name + "_" + station_id
#
#
# ## Step 5: Converts missing data to standard format (this needs to come first because of CWOP_SR missing data flags)
# all_vars = [sfcWind_dir_raw, sfcWind_raw, sfcWind_gust_raw, tas_raw, pr_raw, P_raw, r_raw, p_raw, b_raw, rsds_raw, hurs_raw]
# for item in all_vars:
#     if item == "..." or "....":     # air pressure is 4 chars
#         item = np.nan
#     else:
#         item = float(item)      # converts string text to float
#
# vars_to_keep = ['lat_raw',
#                 'lon_raw',
#                 'sfcWind_raw',  # wind speed
#                 'sfcWind_dir_raw', # wind direction
#                 'tas_raw',  # air temperature
#                 'pr_raw', #
#                 'r_raw', # precipitation -- accumulation in 1 hour
#                 'b_raw', # station pressure
#                 'rsds_raw', # solar radiation
#                 'hurs_raw'] # relative humidity
#
# vars_to_calculate = ['tdps_raw', 'elev']
#
# # Remove these variables from the list of all variables to get drop list
# dropvars = np.setdiff1d(all_vars, vars_to_keep)
#
#
#
# ## Step 6: Tracks existing qa/qc flag for review
# flag1 = 998 # long sequences of 998 is a solar radiation issue, how long do we set as a duration threshold to flag?
# flag2 = '...' # original flag is "..."
#
#
#
#
# ## Step 7: Merge files by station, and outputs cleaned variables as a single .nc file for an individual network
