"""
This script performs data cleaning for the ASOS/AWOS network for
ingestion into the Historical Observations Platform.
Approach:
(1) Read through variables, and calculates derived priority variables if not observed
(2) Drops unnecessary variables
(3) Converts station metadata to standard format, with unique identifier
(4) Converts data metadata to standard format, and converts units into standard units if not provided in standard units.
(5) Converts missing data to standard format
(6) Tracks existing qa/qc flag for review
(7) Merge files by station, and outputs cleaned variables as a single .nc file for an individual network.
Inputs: Raw data for MARITIME stations, with each file representing a month of a year.
Outputs: Cleaned data for an individual network, priority variables, all times. Organized by station as .nc file.
Reference: https://www.ncei.noaa.gov/data/global-hourly/doc/isd-format-document.pdf
"""

# Step 0: Environment set-up
# Import libraries
import os
import xarray as xr
from datetime import datetime
import re
import numpy as np
import csv
import gzip
import pandas as pd
import math

# Import calc_clean.py.
try:
    import calc_clean
except:
    print("Error importing calc_clean.py")
    pass

# Set envr variables and calculate any needed variables

# Experimental/temporary pre-AWS: Set home directory to be where the .git folder is
# Temporary workaround for setting homedir to be head of git repository.
# Only works when current WD is inside of git repository, note!
if os.path.exists('.git'):
    homedir = os.getcwd()
else:
    os.chdir("..")
    if os.path.exists('.git'):
        homedir = os.getcwd()
    else:
        os.chdir("..")
        if os.path.exists('.git'):
            homedir = os.getcwd()
        else:
            os.chdir("..")
            print(os.getcwd())
            if os.path.exists('.git'):
                homedir = os.getcwd()

#homedir = os.getcwd() # 
workdir = "test_platform/data/1_raw_wx/ASOSAWOS/"
savedir = "test_platform/data/2_clean_wx/ASOSAWOS/"
wecc_terr = 'test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp'
wecc_mar = 'test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp' 

# Set up directory to save files, if it doesn't already exist.
try:
    os.mkdir(savedir) # Make the directory to save data in. Except used to pass through code if folder already exists.
except:
    pass

## FUNCTION: Clean ASOS and AWOS data.
# Input: 
# homedir: path to git repository.
# workdir: path to where raw data is saved as .gz files, with each file representing a station's records for 1 year.
# savedir: path to where cleaned files should be saved.
# Output: Cleaned data for an individual network, priority variables, all times. Organized by station as .nc file.

def clean_asosawos(homedir, workdir, savedir):

    os.chdir(workdir)
    
    files = os.listdir() # Gets list of files in directory to work with
    files = list(filter(lambda f: f.endswith(".gz"), files)) # Get list of file names
    
    # Set up error handling.
    errors = {'File':[], 'Time':[], 'Error':[]} # Set up error handling.
    end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download: for error handling csv.
    timestamp = datetime.now().strftime("%m-%d-%Y, %H:%M:%S") # For attributes of netCDF file.

    # Set up lat/lon bounds for filtering data
    try:
        wecc_terr = 'test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp' # Harded-coded because these should not move.
        wecc_mar = 'test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp' 
        t, m, bbox = calc_clean.get_wecc_poly(wecc_terr, wecc_mar)
        lonmin, lonmax = float(bbox['minx']), float(bbox['maxx']) 
        latmin, latmax = float(bbox['miny']), float(bbox['maxy']) 
    except: # If geospatial call fails, hardcode.
        lonmin, lonmax = -139.047795, -102.03721
        latmin, latmax = 30.142739, 60.003861


    ids = list()
    for file in files:
        id = file[0:12] # Remove all YYYYMM (and trailing metadata/file extension) to get file ID.
        if id not in ids:
            ids.append(id)


    # Use ID to grab all files linked to station.
    #for id in ids:
    for id in ['720619-99999', '722689-99999', '726664-94173']: # For testing, pick 3 stations.
        subfiles = list(filter(lambda f: id in f, files))
        subfiles = sorted(subfiles) # Sort files by year in order to concatenate in correct order.
        print(subfiles) # For testing: to know which files are getting 
        
        if bool(subfiles) is False: # If ID returns no files, go to next ID.
            continue

        station_id = "ASOSAWOS_"+id.replace("-", "")
        #print(id, subfiles)
        
        # Initialize dataframe.
        df = pd.DataFrame(columns = ['station_id', 'time', 'latitude', 'longitude', 'elevation', 'qaqc_process', 'ps', 'ps_qc', 'ps_altimeter', 'ps_altimeter_qc', 'psl', 'psl_qc', 'tas', 'tas_qc', 'tdps', 'tdps_qc', 'pr', 'pr_qc', 'pr_duration', 'pr_depth_qc', 'hurs', 'hurs_qc', 'hurs_flag', 'hurs_duration', 'hurs_temp', 'hurs_temp_qc', 'hurs_temp_flag', 'rsds', 'rsds_duration', 'rsds_qc', 'rsds_flag', 'sfcWind', 'sfcWind_qc', 'sfcWind_dir', 'sfcWind_method', 'sfcWind_dir_qc'])
        
        for file in subfiles: # For each file
            #print(file, df)
            try:
                with gzip.open(file, "rt") as f: # Open file
                    csvreader = csv.reader(f)
                    for row in csvreader: # Each row is a record
                        
                        # Initialize all variables and set to be NA by default.
                        string = row[0] # Unpack list.

                        # Filter station by latitude and longitude. Only keep obvs in WECC.
                        # POS 29-34: GEOPHYSICAL-POINT-OBSERVATION latitude coordinate
                        latitude = float(string[28:34])/1000 # Unit: degree
                        # POS 35-41: GEOPHYSICAL-POINT-OBSERVATION longitude coordinate
                        longitude = float(string[34:41])/1000 # Unit: degree
                        if latitude > latmax or latitude < latmin or longitude > lonmax or longitude < lonmin:
                            print("Station {} not in WECC".format([latitude,longitude]))
                            errors['File'].append(file)
                            errors['Time'].append(end_api)
                            errors['Error'].append("File not in WECC. Lat: {} Lon: {}".format(latitude, longitude))
                            break # Go to next file.
                        else:
                            # POS 16-23: GEOPHYSICAL-POINT-OBSERVATION date, # POS 24-27: GEOPHYSICAL-POINT-OBSERVATION time
                            time = datetime.strptime(string[15:27], '%Y%m%d%H%M')
                            
                            # POS 47-51: GEOPHYSICAL-POINT-OBSERVATION elevation dimension
                            elevation = float(string[46:51])
                            # POS 57-60: METEOROLOGICAL-POINT-OBSERVATION quality control process name
                            qaqc_process = string[56:60] 
                            #V01 = No A or M Quality Control applied
                            #V02 = Automated Quality Control
                            #V03 = subjected to Quality Control

                            # Mandatory data

                            #POS: 61-63: WIND-OBSERVATION direction angle
                            sfcWind_dir = int(string[60:63]) # Units degrees

                            #POS: 64-64: WIND-OBSERVATION direction quality code
                            sfcWind_dir_qc = string[63]

                            #POS: 65-65 WIND-OBSERVATION type code
                            sfcWind_method = string[64]

                            #POS: 66-69: WIND-OBSERVATION speed rate
                            sfcWind = float(string[65:69])/10 # Units m/s

                            #POS: 70-70: WIND-OBSERVATION speed quality code
                            sfcWind_qc = string[69]
                            # Note: One row returns A here, this is an error.

                            #POS 88-92: AIR-TEMPERATURE-OBSERVATION air temperature
                            tas = float(string[87:92])/10

                            #POS 93: AIR-TEMPERATURE-OBSERVATION air temperature quality code
                            tas_qc = string[92]

                            #POS 94-98: AIR-TEMPERATURE-OBSERVATION dew point temperature
                            tdps = float(string[93:98])/10

                            #POS 99-99: AIR-TEMPERATURE-OBSERVATION dew point quality code
                            tdps_qc = string[98]

                            #POS 100-104: ATMOSPHERIC-PRESSURE-OBSERVATION sea level pressure
                            psl = float(string[99:104])/10 # In hectopascals, CONVERT.

                            #POS 105-105: ATMOSPHERIC-PRESSURE-OBSERVATION sea level pressure quality code
                            psl_qc = string[104]
                            
                            # Additional data - begins with ADD
                            
                            # Liquid precipitation - begins with AA[1-4]

                            # Figure out length of precip. string.
                            precip_len = re.search("(?<=AA)[\d]", string)

                            if precip_len is not None: # If precip data exists
                                # Get everything after AA#
                                # QA/QC codes include letters. So, use precip_length to pull string of correct length.
                                precip = re.search("(?<=AA1|AA2|AA3|AA4)[\da-zA-Z]{8}", string).group() # Get the 8 strings after AA1/2/3/4
                                pr_duration = int(precip[0:2]) # Hours
                                pr = float(precip[2:6])/10 # In mm
                                pr_depth_qc = int(precip[6:7]) # QA for precipitation depth
                                pr_qc = precip[7:8] # QA for precipitation data.
                                
                                # Take first measurement as primary, unless missing.
                                if float(precip[2:6]) == 9999:
                                    if str(precip_len.group()) != "1":
                                        # If precip depth of first report is missing.
                                        precip = re.search("(?<=AA1|AA2|AA3|AA4)[\da-zA-Z]{16}", string).group() # Get second set of precip values.
                                        pr_duration = int(precip[8:10]) # Hours
                                        pr = float(precip[10:14])/10 # In mm
                                        pr_depth_qc = int(precip[14:15]) # QA for precipitation depth
                                        pr_qc = precip[15:16] # QA for precipitation data.
                            else:
                                pr = np.nan
                                pr_qc = np.nan
                                pr_duration = np.nan
                                pr_depth_qc = np.nan

                            # Relative humidity - shouldn't be in this dataset, but capture if it is.
                            
                            hurs_string = re.search("(?<=CH1|CH2)[\da-zA-Z]{15}", string) #Section starts with CH 
                            
                            if hurs_string is not None: # If precip data exists
                                hurs_string = hurs_string.group() # Access string from match.
                                
                                hurs_duration = int(hurs_string[0:2]) # Minutes
                                hurs_temp = int(hurs_string[2:7])/10 # In deg. C
                                hurs_temp_qc = hurs_string[7]
                                hurs_temp_flag = int(hurs_string[8])
                                hurs = int(hurs_string[9:13])/10 # In percent
                                hurs_qc = hurs_string[13]
                                hurs_flag = int(hurs_string[14])

                            else:
                                hurs = np.nan
                                hurs_qc = np.nan
                                hurs_flag = np.nan
                                hurs_duration = np.nan
                                hurs_temp = np.nan
                                hurs_temp_qc = np.nan
                                hurs_temp_flag = np.nan
                            
                            # Solar radiation - "global irradiance" in our dataset.
                            # Not keeping "direct beam irradiance", "diffuse irradiance" which are also provided in same string.
                            rsds_string = re.search("(?<=GM1)[\da-zA-Z]{11}", string) #Section starts with CH 
                            
                            if rsds_string is not None: # If precip data exists
                                rsds_string = rsds_string.group() # Access string from match.
                                
                                rsds_duration = float(rsds_string[0:4]) # Time period over which solar radiation integrated, in minutes
                                rsds = float(rsds_string[4:8]) # In w/m2
                                rsds_flag = rsds_string[8:10]
                                rsds_qc = rsds_string[10:12]

                            else:
                                rsds = np.nan
                                rsds_duration = np.nan
                                rsds_qc = np.nan
                                rsds_flag = np.nan

                            # Station pressure
                            ps_string = re.search("(?<=MA1)[\da-zA-Z]{12}", string) #Section starts with CH 
                            
                            if ps_string is not None: # If precip data exists
                                ps_string = ps_string.group() # Access string from match.
                                ps_altimeter = float(ps_string[0:5])/10 # HPa
                                ps_altimeter_qc = ps_string[5]
                                ps = float(ps_string[6:11])/10 # HPa
                                ps_qc = ps_string[11]

                            else:
                                ps = np.nan
                                ps_qc = np.nan
                                ps_altimeter = np.nan
                                ps_altimeter_qc = np.nan

                            # Standardize NAs
                            try:
                                if latitude == 99.999 or longitude == 99.999:
                                    continue # If lat or lon values are NA, skip observation.
                                if elevation == 9999:
                                    elevation = np.nan
                                if ps == 9999.9:
                                    ps = np.nan
                                if ps_altimeter == 9999.9:
                                    ps_altimeter = np.nan
                                if psl == 9999.9:
                                    psl = np.nan
                                if tas == 999.9:
                                    tas = np.nan
                                if tdps == 999.9:
                                    tdps = np.nan
                                if pr == 999.9:
                                    pr = np.nan
                                if pr_duration == 99:
                                    pr_duration = np.nan
                                if hurs == 999.9:
                                    hurs = np.nan
                                if hurs_duration == 99:
                                    hurs_duration = np.nan
                                if hurs_temp == 999.9:
                                    hurs_temp = np.nan
                                if rsds == 9999:
                                    rsds = np.nan
                                if rsds_duration == 9999:
                                    rsds_duration = np.nan
                                if sfcWind_dir == 999:
                                    sfcWind_dir = np.nan
                                if sfcWind == 999.9:
                                    sfcWind = np.nan
                            except Exception as e:
                                print(e) # Could add more robust error handling here if desired, but should not be necessary.
                            
                            # For each row of data, append data to row.
                            data = [station_id, time, latitude, longitude, elevation, qaqc_process, ps, ps_qc, ps_altimeter, ps_altimeter_qc, psl, psl_qc, tas, tas_qc, tdps, tdps_qc, pr, pr_qc, pr_duration, pr_depth_qc, hurs, hurs_qc, hurs_flag, hurs_duration, hurs_temp, hurs_temp_qc, hurs_temp_flag, rsds, rsds_duration, rsds_qc, rsds_flag, sfcWind, sfcWind_qc, sfcWind_dir, sfcWind_method, sfcWind_dir_qc]
                            data = pd.DataFrame([data], columns = df.columns)
                            df = df.append(data, ignore_index = True)
                            
                            # For testing: progress update. Print status update every 1k rows.
                            # If station reports every 10 min, expect up to 50k observations per year.
                            if len(df)%1000==0:
                                print("{} observations parsed.".format(len(df)))
                            
            except Exception as e:
                print(file, e)
                errors['File'].append(file)
                errors['Time'].append(end_api)
                errors['Error'].append(e)

        if df.empty is True:
            print("df {} not saved".format(file))

        if df.empty is False: # If there is data in the dataframe, convert to xarray object.
            try:
                ds = df.to_xarray()
                # Update dimensions and coordinates

                # Add dimensions: station ID and time.
                ds = ds.set_coords('time').swap_dims({'index': 'time'}) # Swap index with time.
                ds = ds.assign_coords(id = str(station_id))
                ds = ds.expand_dims("id") # Add station_id as index.
                ds = ds.drop_vars(("station_id", "index")) # Drop station_id variable and index coordinate.
                ds = ds.rename({'id': 'station'}) # Rename id to station_id.
                
                # Add coordinates: latitude and longitude.
                ds = ds.set_coords(("latitude", "longitude"))

                # Update dimension and coordinate attributes.

                # Time
                ds['time'].attrs['long_name'] = "time"
                ds['time'].attrs['standard_name'] = "time"
                ds['time'].attrs['comment'] = "In UTC."
                
                # Station ID
                ds['station'].attrs['long_name'] = "station_id"
                ds['station'].attrs['comment'] = "Unique ID created by Eagle Rock Analytics. Includes network name appended to original unique station ID provided by network."
                
                # Latitude
                ds['latitude'].attrs['long_name'] = "latitude"
                ds['latitude'].attrs['standard_name'] = "latitude"
                ds['latitude'].attrs['units'] = "degrees_north"
                
                # Longitude
                ds['longitude'].attrs['long_name'] = "longitude"
                ds['longitude'].attrs['standard_name'] = "longitude"
                ds['longitude'].attrs['units'] = "degrees_east"

                # Elevation
                ds['elevation'].attrs['long_name'] = "station_elevation"
                ds['elevation'].attrs['units'] = "meter"
                ds['elevation'].attrs['positive'] = "up" # Define which direction is positive
                
                # Update global attributes
                ds = ds.assign_attrs(title = "ASOS/AWOS cleaned")
                ds = ds.assign_attrs(institution = "Eagle Rock Analytics / Cal Adapt")
                ds = ds.assign_attrs(source = "")
                ds = ds.assign_attrs(history = "ASOSAWOS_clean.py script run on {} UTC".format(timestamp))
                ds = ds.assign_attrs(comment = "Intermediate data product: may not have been subject to any cleaning or QA/QC processing")
                ds = ds.assign_attrs(license = "")
                ds = ds.assign_attrs(citation = "")
                ds = ds.assign_attrs(disclaimer = "This document was prepared as a result of work sponsored by the California Energy Commission (PIR-19-006). It does not necessarily represent the views of the Energy Commission, its employees, or the State of California. Neither the Commission, the State of California, nor the Commission's employees, contractors, or subcontractors makes any warranty, express or implied, or assumes any legal liability for the information in this document; nor does any party represent that the use of this information will not infringe upon privately owned rights. This document has not been approved or disapproved by the Commission, nor has the Commission passed upon the accuracy of the information in this document.")
                
                # Update variable attributes and do unit conversions

                # tas: air surface temperature (K)
                if "tas" in ds.keys():
                    ds['tas_raw'] = ds['tas'] # Move original data to raw column.
                    try: 
                        ds['tas'] = calc_clean._unit_degC_to_K(ds['tas'])
                    except:
                        print("tas: calc_clean.py not working.")
                        ds['tas'] = ds['tas'] + 273.15 # Convert to K (backup method)
                    ds['tas'].attrs['long_name'] = "air_temperature"
                    ds['tas'].attrs['standard_name'] = "air_temperature"
                    ds['tas'].attrs['units'] = "degree_Kelvin"
                    ds['tas'].attrs['ancillary_variables'] = "tas_raw tas_qc" # List other variables associated with variable (QA/QC)
                    ds['tas'].attrs['comment'] = "Converted from Celsius."
                    
                    ds['tas_raw'].attrs['long_name'] = "air_temperature"
                    ds['tas_raw'].attrs['standard_name'] = "air_temperature"
                    ds['tas_raw'].attrs['units'] = "degree_Celsius"
                    ds['tas_raw'].attrs['ancillary_variables'] = "tas tas_qc" # List other variables associated with variable (QA/QC)

                if "tas_qc" in ds.keys():
                    ds['tas_qc'].attrs['flag_values'] = "[0, 1, 2, 3, 4, 5, 6, 7, 9, A, C, I, M, P, R, U]"
                    ds['tas_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                    
                # ps: surface air pressure (Pa)
                if 'ps' in ds.keys():
                    ds['ps_raw'] = ds['ps'] # Move original data to raw column.
                    ds['ps'].attrs['long_name'] = "station_air_pressure"
                    ds['ps'].attrs['standard_name'] = "air_pressure"
                    ds['ps'].attrs['units'] = "hPa"
                    ds['ps'].attrs['ancillary_variables'] = "ps_raw ps_qc" # List other variables associated with variable (QA/QC)
                    
                    ds['ps_raw'].attrs['long_name'] = "station_air_pressure"
                    ds['ps_raw'].attrs['standard_name'] = "air_pressure"
                    ds['ps_raw'].attrs['units'] = "hPa"
                    ds['ps_raw'].attrs['ancillary_variables'] = "ps ps_qc ps_altimeter ps_altimeter_qc" # List other variables associated with variable (QA/QC)

                    # Delete sea level pressure if station pressure included.
                    ds.drop(("psl", "psl_qc"))

                if "ps_qc" in ds.keys():
                    ds['ps_qc'].attrs['flag_values'] = "[0, 1, 2, 3, 4, 5, 6, 7, M, 9]"
                    ds['ps_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
        
                if "ps_altimeter" in ds.keys():
                    ds['ps_altimeter'].attrs['long_name'] = "altimeter_setting"
                    ds['ps_altimeter'].attrs['units'] = "hPa"
                    ds['ps_altimeter'].attrs['ancillary_variables'] = "ps ps_qc ps_altimeter ps_altimeter_qc" # List other variables associated with variable (QA/QC)
                    ds['ps_altimeter'].attrs['comment'] = "The pressure value to which an aircraft altimeter is set so that it will indicate the altitude relative to mean sea level of an aircraft on the ground at the location for which the value was determined." # Description of variable meaning, as not CF-standard.
                
                if "ps_altimeter_qc" in ds.keys():
                    ds['ps_altimeter_qc'].attrs['flag_values'] = "[0, 1, 2, 3, 4, 5, 6, 7, M, 9]"
                    ds['ps_altimeter_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."

                # If station air pressure not reported, convert from sea level pressure.
                # Inputs: sea level pressure (mb/hPa), elevation (m), and air temperature (K)
                elif 'psl' in ds.keys():
                    ds['psl_raw'] = ds['ps'] # Move original data to raw column.
                    
                    # Convert from sea level pressure. TO DO. (Need to know elev of barometer.)
                    try:
                        ds['ps'] = calc_clean._calc_ps(ds['psl'], ds['elevation'], ds['tas'])
                    except:
                        print("ps: calc_clean.py not working.") 
                        ds['ps'] = ds['psl']*math.e**(-ds['elevation']/(ds['tas']*29.263)) # Backup method, remove once calc_clean works consistently?

                    ds['ps'].attrs['long_name'] = "station_air_pressure"
                    ds['ps'].attrs['standard_name'] = "air_pressure"
                    ds['ps'].attrs['units'] = "hPa"
                    ds['ps'].attrs['ancillary_variables'] = "psl_raw psl_qc elevation" # List other variables associated with variable (QA/QC)
                    ds['ps'].attrs['comment'] = "Converted from sea level pressure."
                    
                    ds['psl_raw'].attrs['long_name'] = "sea_level_air_pressure"
                    ds['psl_raw'].attrs['standard_name'] = "air_pressure_at_sea_level"
                    ds['psl_raw'].attrs['units'] = "hPa"
                    ds['psl_raw'].attrs['ancillary_variables'] = "ps psl_qc" # List other variables associated with variable (QA/QC)
                
                # tdps: dew point temperature (K)
                # tdps always provided by ISD reports, so no need to calculate here.
                if "tdps" in ds.keys(): # If variable already exists, rename.
                    ds['tdps_raw'] = ds['tdps'] # Move original data to raw column.

                    try:
                        ds['tdps'] = calc_clean._unit_degC_to_K(ds['tdps'])
                    except:
                        print("tdps: calc_clean.py not working.") 
                        ds['tdps'] = ds["tdps"]+273.15 # Convert from C to K
                    ds['tdps'].attrs['long_name'] = "dew_point_temperature"
                    ds['tdps'].attrs['standard_name'] = "dew_point_temperature"
                    ds['tdps'].attrs['units'] = "degree_Kelvin"
                    ds['tdps'].attrs['ancillary_variables'] = "tdps_raw tdps_qc" # List other variables associated with variable (QA/QC)
                    ds['tdps'].attrs['comment'] = "Converted from Celsius."
                    
                    ds['tdps_raw'].attrs['long_name'] = "dew_point_temperature"
                    ds['tdps_raw'].attrs['standard_name'] = "dew_point_temperature"
                    ds['tdps_raw'].attrs['units'] = "degree_Celsius"
                    ds['tdps_raw'].attrs['ancillary_variables'] = "tdps tdps_qc" # List other variables associated with variable (QA/QC)

                if 'tdps_qc' in ds.keys():
                    ds['tdps_qc'].attrs['flag_values'] = "[0, 1, 2, 3, 4, 5, 6, 7, 9, A, C, I, M, P, R, U]"
                    ds['tdps_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                    
                # pr: precipitation (Flag for discussion about CF compliance and units) 
                if 'pr' in ds.keys():
                    ds['pr_raw'] = ds['pr'] # Move original data to raw column.
                    # Note leave final conversion here to next stage.
                    ds['pr'].attrs['long_name'] = "precipitation"
                    ds['pr'].attrs['standard_name'] = "" # Flag for discussion.
                    ds['pr'].attrs['units'] = "" # Flag for discussion.
                    ds['pr'].attrs['ancillary_variables'] = "pr_raw pr_qc pr_depth_qc pr_duration" # List other variables associated with variable (QA/QC)
                    ds['pr'].attrs['comment'] = "" # To be completed.

                    ds['pr_raw'].attrs['long_name'] = "precipitation"
                    ds['pr_raw'].attrs['standard_name'] = "" # TBD. Maybe delete, not standard.
                    ds['pr_raw'].attrs['units'] = "mm"
                    ds['pr_raw'].attrs['ancillary_variables'] = "pr pr_qc pr_depth_qc pr_duration" # List other variables associated with variable (QA/QC)

                    ds['pr_duration'].attrs['long_name'] = "precipitation measurement interval"
                    ds['pr_duration'].attrs['units'] = "hours"
                    ds['pr_duration'].attrs['ancillary_variables'] = "pr pr_raw pr_qc pr_depth_qc" # List other variables associated with variable (QA/QC)

                if 'pr_qc' in ds.keys():
                    ds['pr_qc'].attrs['flag_values'] = "[0, 1, 2, 3, 4, 5, 6, 7, 9, A, I, M, P, R, U]"
                    ds['pr_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                    ds['pr_depth_qc'].attrs['flag_values'] = "[1, 2, 3, 4, 5, 6, 7, 8, E, I, J, 9]"
                    ds['pr_depth_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                
                # hurs: relative humidity
                if 'hurs' in ds.keys():
                    ds['hurs'].attrs['long_name'] = "average_relative_humidity"
                    ds['hurs'].attrs['standard_name'] = "relative_humidity" 
                    ds['hurs'].attrs['units'] = "percent" 
                    ds['hurs'].attrs['ancillary_variables'] = "hurs_qc hurs_flag hurs_duration hurs_temp hurs_temp_qc hurs_temp_flag" # List other variables associated with variable (QA/QC)

                    ds['hurs_duration'].attrs['long_name'] = "relative_humidity_duration"
                    ds['hurs_duration'].attrs['units'] = "minutes" 
                    ds['hurs_duration'].attrs['ancillary_variables'] = "hurs hurs_qc hurs_flag hurs_temp hurs_temp_qc hurs_temp_flag" # List other variables associated with variable (QA/QC)

                    ds['hurs_temp'].attrs['long_name'] = "average_air_temperature_at_hurs_instrument"
                    ds['hurs_temp'].attrs['units'] = "degree_celsius" 
                    ds['hurs_temp'].attrs['ancillary_variables'] = "hurs hurs_qc hurs_flag hurs_temp hurs_temp_qc hurs_temp_flag" # List other variables associated with variable (QA/QC)
                
                elif 'tas' and 'tdps' in ds.keys():
                    try:
                        ds['hurs'] = calc_clean._calc_relhumid(ds['tas'], ds['tdps'])
                    except:
                        print("hurs: calc_clean.py not working.")
                        es = 0.611 * np.exp(5423 * ((1/273) - (1/tas)))
                        e = 0.611 * np.exp(5423 * ((1/273) - (1/tdps)))
                        ds['hurs'] = 100 * (e/es)
                    ds['hurs'].attrs['long_name'] = "relative_humidity" # Note some here will be avg. and some not, will have to align when we merge down.
                    ds['hurs'].attrs['standard_name'] = "relative_humidity" 
                    ds['hurs'].attrs['units'] = "percent" 
                    ds['hurs'].attrs['comment'] = "Calculated from air temperature and dew point temperature."

                if 'hurs_qc' in ds.keys():
                    ds['hurs_qc'].attrs['flag_values'] = "[1, 3, 9]"
                    ds['hurs_qc'].attrs['flag_meanings'] =  "passed_all_qc_checks failed_all_qc_checks missing"
                    
                    ds['hurs_flag'].attrs['flag_values'] = "[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]"
                    ds['hurs_flag'].attrs['flag_meanings'] =  "See QA/QC csv for network."

                    ds['hurs_temp_qc'].attrs['flag_values'] = "[1, 3, 9]"
                    ds['hurs_temp_qc'].attrs['flag_meanings'] =  "passed_all_qc_checks failed_all_qc_checks missing"
                    
                    ds['hurs_temp_flag'].attrs['flag_values'] = "[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]"
                    ds['hurs_temp_flag'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                    
                # rsds: surface_downwelling_shortwave_flux_in_air (solar radiation, w/m2)
                if 'rsds' in ds.keys():
                    #print(ds['rsds'].values[0])
                    if np.isnan(ds['rsds'].values[0]).all(): # If all values are NA
                        pass # Continue and don't add this variable.
                    else: # Otherwise, add variables.
                        ds['rsds'].attrs['long_name'] = "solar_radiation"
                        ds['rsds'].attrs['standard_name'] = "surface_downwelling_shortwave_flux_in_air"
                        ds['rsds'].attrs['units'] = "W m-2"
                        ds['rsds'].attrs['ancillary_variables'] = "rsds_duration rsds_qc rsds_flag" # List other variables associated with variable (QA/QC)
                        ds['rsds'].attrs['comment'] = "Waveband ranges from 0.4 - 2.3 micrometers."

                        ds['rsds_duration'].attrs['long_name'] = "time_period"
                        ds['rsds_duration'].attrs['units'] = "minutes"
                        ds['rsds_duration'].attrs['ancillary_variables'] = "rsds rsds_qc rsds_flag" # List other variables associated with variable (QA/QC)
                        ds['rsds_duration'].attrs['comment'] = "Time period over which solar radiation is integrated."

                    if 'rsds_qc' in ds.keys():
                        ds['rsds_qc'].attrs['flag_values'] = "[0, 1, 2, 3, 9]"
                        ds['rsds_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                    
                        ds['rsds_flag'].attrs['flag_values'] = "00, 01, 02, 03, 04, 05, 06, 07, 08, 09, 10-93, 94-97, 98, 99"
                        ds['rsds_flag'].attrs['flag_meanings'] =  "See QA/QC csv for network."

                   
                # sfcWind : wind speed (m/s) (Method of calculation may vary, standardize during hourly merge or QA/QC process.)
                if "sfcWind" in ds.keys(): # No conversions needed, do not add raw column.
                    ds['sfcWind'].attrs['long_name'] = "wind_speed"
                    ds['sfcWind'].attrs['standard_name'] = "wind_speed"
                    ds['sfcWind'].attrs['units'] = "m s-1"
                    ds['sfcWind'].attrs['ancillary_variables'] = "sfcWind_qc sfcWind_method" # List other variables associated with variable (QA/QC)
                    ds['sfcWind'].attrs['comment'] = "Method of wind speed calculation varies, see sfcWind_method."

                if 'sfcWind_method' in ds.keys():
                    ds['sfcWind_method'].attrs['long_name'] = "wind_speed_calculation_method"
                    ds['sfcWind_method'].attrs['flag_values'] = "[A, B, C, H, N, R, Q, T, V, 9]"
                    ds['sfcWind_method'].attrs['flag_meanings'] = "abridged_beaufort beaufort calm 5-minute_average normal 60-minute_average squall 180-minute_average variable missing" 

                if 'sfcWind_qc' in ds.keys():
                    ds['sfcWind_method'].attrs['flag_values'] = "[0, 1, 2, 3, 4, 5, 6, 7, 9]"
                    ds['sfcWind_method'].attrs['flag_meanings'] = "See QA/QC csv for network." 


                # sfcWind_dir: wind direction
                if "sfcWind_dir" in ds.keys(): # No conversions needed, do not make raw column.
                    ds['sfcWind_dir'].attrs['long_name'] = "wind_direction"
                    ds['sfcWind_dir'].attrs['standard_name'] = "wind_from_direction"
                    ds['sfcWind_dir'].attrs['units'] = "degrees_clockwise_from_north"
                    ds['sfcWind_dir'].attrs['ancillary_variables'] = "sfcWind_dir_qc" # List other variables associated with variable (QA/QC)
                    
                if 'sfcWind_dir_qc' in ds.keys():
                    ds['sfcWind_dir_qc'].attrs['flag_values'] = "[0, 1, 2, 3, 4, 5, 6, 7, 9]"
                    ds['sfcWind_dir_qc'].attrs['flag_meanings'] = "See QA/QC csv for network." 


                # Update attributes for any non-standard variables
                
                #QAQC process
                if "qaqc_process" in ds.keys():
                    ds['qaqc_process'].attrs['long_name'] = "qaqc_process_type"
                    ds['qaqc_process'].attrs['flag_values'] = "[V01, V02, V03]"
                    ds['qaqc_process'].attrs['flag_meanings'] = "no_qaqc automated_qaqc subjected_to_qaqc" 

                # Data source
                if "data_source" in ds.keys():
                    ds['data_source'].attrs['long_name'] = "source_of_data"
                    ds['data_source'].attrs['flag_values'] = "[6, 7]"
                    ds['data_source'].attrs['flag_meanings'] = "ASOSAWOS ASOSAWOS_merged_with_USAF_surface_hourly" 

                # # Reorder variables
                # In following order:
                desired_order = ['ps', 'tas', 'tdps', 'pr', 'hurs', 'rsds', 'sfcWind', 'sfcWind_dir']
                desired_order = [i for i in desired_order if i in list(ds.keys())] # Only keep vars which are in ds.
                rest_of_vars = [i for i in list(ds.keys()) if i not in desired_order] # Retain rest of variables at the bottom.
                new_index = desired_order + rest_of_vars
                ds = ds[new_index]

                # Testing: Manually check values to see that they seem correctly scaled, no unexpected NAs.
                # for var in ds.variables:
                #     try:
                #         print([var, float(ds[var].min()), float(ds[var].max())]) 
                #     except:
                #         next
                
            except Exception as e: # If error in xarray reorganization
                print(file, e)
                errors['File'].append(station_id)
                errors['Time'].append(end_api)
                errors['Error'].append(e)
                    

            # Write station file to netcdf.
            if ds is None: # Should be caught be error handling above, but add in case.
                print("ds {} not saved.".format(file))
                continue
            
            else:
                try:
                    print(ds) # For testing.
                    filename = station_id+".nc" # Make file name
                    filepath = homedir+"/"+savedir+filename # Write file path
                    print(filepath)

                    ds.to_netcdf(path = filepath) # Save station file.
                    print("Saving {} with dims {}".format(filename, ds.dims))
                    ds.close() # Close dataframe.
                except Exception as e:
                    print(filename, e)
                    errors['File'].append(filename)
                    errors['Time'].append(end_api)
                    errors['Error'].append(e)
                    continue            

    #Write errors to csv
    filepath = homedir+"/"+savedir+"errors_asosawos_{}.csv".format(end_api) # Set path to save error file.
    print(filepath)
    #print(errors)
    with open(filepath, "w") as outfile:
        # pass the csv file to csv.writer function.
        writer = csv.writer(outfile)

        # pass the dictionary keys to writerow
        # function to frame the columns of the csv file
        writer.writerow(errors.keys())

        # make use of writerows function to append
        # the remaining values to the corresponding
        # columns using zip function.
        writer.writerows(zip(*errors.values()))

# Run function
clean_asosawos(homedir, workdir, savedir)

# Testing:
## Import file.
# os.chdir(savedir)
# test = xr.open_dataset("ASOSAWOS_72061999999.nc")
# print(test) 

# # # # Test 1: multi-year merges work as expected.
# print(str(test['time'].min())) # Get start time
# print(str(test['time'].max())) # Get end time


# # ## Test 2: Inspect vars and attributes
# # ## 
# for var in test.variables: 
#     try:
#         print([var, float(test[var].min()), float(test[var].max())]) 
#     except:
#         continue

# # # Test 3: Get one month's data and test subsetting.
# print(test.sel(time = "2015-05"))

# # # Next file.
# test = xr.open_dataset("ASOSAWOS_72268999999.nc") 
# print(test)
# print(str(test['time'].min())) # Get start time
# print(str(test['time'].max())) # Get end time

# ## Inspect vars and attributes
# ## 
# for var in test.variables: 
#     try:
#         print([var, float(test[var].min()), float(test[var].max())]) 
#     except:
#         continue

# # # Get a few rows
# print(test.sel(time = "1989-06")) # Should only be 3 observations here. Very incomplete year.
