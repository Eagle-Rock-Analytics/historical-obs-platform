"""
This is a script where Stage 3: QA/QC related common functions, conversions, and operations is stored for ease of use
for the Historical Observations Platform.
"""

## Import Libraries
import boto3
import geopandas as gp
import numpy as np
import pandas as pd
import requests
import urllib
import datetime
import math


## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client('s3') # for lower-level processes

## Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"


#----------------------------------------------------------------------
## QA/QC Helper functions
def get_file_paths(network):
    rawdir = "1_raw_wx/{}/".format(network)
    cleandir = "2_clean_wx/{}/".format(network)
    qaqcdir = "3_qaqc_wx/{}/".format(network)
    mergedir = "4_merge_wx/{}/".format(network)
    return rawdir, cleandir, qaqcdir, mergedir


def get_wecc_poly(terrpath, marpath):
    """
    Identifies a bbox of WECC area to filter stations against
    Input vars: shapefiles for maritime and terrestrial WECC boundaries
    Returns: spatial objects for each shapefile, and bounding box for their union.
    """
    t = gp.read_file(terrpath)  ## Read in terrestrial WECC shapefile.
    m = gp.read_file(marpath)   ## Read in marine WECC shapefile.
    bbox = t.union(m).bounds    ## Combine polygons and get bounding box of union.
    return t,m, bbox


#----------------------------------------------------------------------
## Part 1 functions (whole station/network)

## missing spatial coords (lat-lon)
def qaqc_missing_latlon(df):
    """
    Checks if latitude and longitude is missing for a station.
    If missing, station is flagged to not proceed through QA/QC.
    """

    # latitude
    if df['lat'].isnull().values.any() == True:
        df = None # returns empty df to flag that it does not pass
    else:
        df = df

    # longitude
    if df['lon'].isnull().values.any() == True:
        df = None # returns empty df to flag that it does not pass
    else:
        df = df

    return df

## in bounds of WECC
def qaqc_within_wecc(df):
    """
    Checks if station is within terrestrial & marine WECC boundaries.
    If outside of boundaries, station is flagged to not proceed through QA/QC.
    """

    t, m, bbox = get_wecc_poly(wecc_terr, wecc_mar) # Call get_wecc_poly
    lat_to_check = df['lat'].iloc[0]
    lon_to_check = df['lon'].iloc[0]

    # latitude
    if (lat_to_check < bbox.miny.values) or (lat_to_check > bbox.maxy.values):
        df = pd.DataFrame() # returns empty df to flag that it does not pass
    else:
        df = df

    # longitude
    if (lon_to_check > bbox.maxx.values) or (lon_to_check < bbox.minx.values):
        df = pd.DataFrame() # returns empty df to flag that it does not pass
    else:
        df = df

    return df

## elevation
def _grab_dem_elev_m(lat_to_check, lon_to_check):
    """
    Pulls elevation value from the USGS Elevation Point Query Service, lat lon must be in decimal degrees (which it is after cleaning)
    Modified from: https://gis.stackexchange.com/questions/338392/getting-elevation-for-multiple-lat-long-coordinates-in-python
    """

    url = r'https://epqs.nationalmap.gov/v1/json?'

    # define rest query params
    params = {
        'output': 'json',
        'x': lon_to_check,
        'y': lat_to_check,
        'units': 'Meters'
    }

    # format query string and return value
    result = requests.get((url + urllib.parse.urlencode(params)))
    dem_elev_long = float(result.json()['value'])
    dem_elev_short = '{:.2f}'.format(dem_elev_long) # make sure to round off lat-lon values so they are not improbably precise for our needs

    return dem_elev_short


def qaqc_elev_infill(df):
    """
    Checks if elevation is NA/missing. If missing, fill in elevation from either DEM or station.
    Some stations have all nan elevation values (e.g., NDBC, MARITIME)
    Some stations have single/few but not all nan elevation values (e.g., otherisd, asosawos)
    """

    print('Elevation values pre-infilling: {}'.format(df['elevation'].unique()))
    print('Elevation eraqc values pre-infilling: {}'.format(df['elevation_eraqc'].unique())) # testing

    # first check to see if any elev value is missing
    if df['elevation'].isnull().any() == True: 

        # all elevation values are reported as nan (some ndbc/maritime)
        if df['elevation'].isnull().values.all() == True: 
            # print('This station reports all missing elevation -- infilling from DEM') # testing

            try:  # in-fill if value is missing
                # check if lat-lon has changed over time
                nan_lats = df['lat'].unique()
                nan_lons = df['lon'].unique()

                if (len(nan_lats) == 1) and (len(nan_lons) == 1): # single lat-lon pair for missing elevs
                    try:
                        dem_elev_value = _grab_dem_elev_m(df['lat'].iloc[0], df['lon'].iloc[0])
                        df.loc[df['elevation'].isnull(), 'elevation_eraqc'] = 3 # see era_qaqc_flag_meanings.csv
                        df.loc[df['elevation'].isnull(), 'elevation'] = float(dem_elev_value)
                    except: # some buoys out of range of dem (past coastal range) report nan elevation, manually set to 0.00m and flag
                        df.loc[df['elevation'].isnull(), 'elevation_eraqc'] = 5 # see era_qaqc_flag_meanings.csv
                        df.loc[df['elevation'].isnull(), 'elevation'] = float(0.00) # manual infilling

                else: # multiple pairs of lat-lon for missing elevs
                    for ilat in nan_lats:
                        for ilon in nan_lons:
                            dem_elev_value = _grab_dem_elev_m(ilat, ilon)
                            df.loc[(df['lat'] == ilat) & (df['lon'] == ilon), 'elevation_eraqc'] = 3 # see era_qaqc_flag_meanings.csv
                            df.loc[(df['lat'] == ilat) & (df['lon'] == ilon), 'elevation'] = float(dem_elev_value)
                        
            except: # elevation cannot be obtained from DEM
                df = pd.DataFrame() # returns empty df to flag that it does not pass


        # some stations have a single/few nan reported (some otherisd/asosawos stations)
        else:   # multiple values for elevation, infill each instance if missing/incorrectly coded (e.g., zeros when shouldnt be)
            # print('This station reports a missing elevation, but not all -- targeted in-filling') # dummy flag message to note which stations this occurs for focused testing

            try:
                # locate all instances of nan values as elevation codes
                nan_coded = df[df['elevation'].isnull()]
                nan_lats = nan_coded['lat'].unique()
                nan_lons = nan_coded['lon'].unique()

                if (len(nan_lats) == 1) and (len(nan_lons) == 1): # single lat-lon pair for missing elevs
                    if (nan_lats[0] == df['lat'].iloc[0]) & (nan_lons[0] == df['lon'].iloc[0]): # single set of lat-lons matches station, infill from station
                        df.loc[df['elevation'].isnull(), 'elevation_eraqc'] = 4 # see era_qaqc_flag_meanings.csv
                        df.loc[df['elevation'].isnull(), 'elevation'] = df['elevation'].iloc[0]

                    else: # lat-lon of missing elev does not match station lat-lon (has shifted), infill from dem
                        dem_elev_value = _grab_dem_elev_m(nan_lats[0], nan_lons[0])
                        df.loc[df['elevation'].isnull(), 'elevation_eraqc'] = 3 # see era_qaqc_flag_meanings.csv
                        df.loc[df['elevation'].isnull(), 'elevation'] = float(dem_elev_value)

                else: # multiple pairs of lat-lon for missing elevs
                    for ilat in nan_lats:
                        for ilon in nan_lons:
                            dem_elev_value = _grab_dem_elev_m(ilat, ilon)
                            df.loc[(df['lat'] == ilat) & (df['lon'] == ilon), 'elevation_eraqc'] = 3 # see era_qaqc_flag_meanings.csv
                            df.loc[(df['lat'] == ilat) & (df['lon'] == ilon), 'elevation'] = float(dem_elev_value)
                                                            
            except: # elevation cannot be obtained from DEM
                df = pd.DataFrame() # returns empty df to flag that it does not pass

    else:
        df = df

    return df


def qaqc_elev_range(df):
    """
    Checks valid values to identify suspicious elevation values that are larger than 10m in difference
    Checks if valid elevation value is outside of range of reasonable values for WECC region.
    If outside range, station is flagged to not proceed through QA/QC.
    """

    # first check for suspicious values
    elev_vals = df['elevation'].unique() # identify how many different elevation "values" are present

    # elevation values flagged as incorrectly coded
    # uses a threshold of 10m different from the station elevation to identify suspicious elevations
    for elev_value in elev_vals:
        if (elev_value > df['elevation'].iloc[0] + 10) or (elev_value < df['elevation'].iloc[0] - 10): # 10m above and below bounds to check
            off_elevs = df.loc[df['elevation'] == elev_value]
            off_lats = off_elevs['lat'].unique()
            off_lons = off_elevs['lon'].unique()

            if (len(off_lats) == 1) and (len(off_lons) == 1): # single lat-lon pair for incorrectly coded elevation
                if (off_lats[0] == df['lat'].iloc[0]) & (off_lons[0] == df['lon'].iloc[0]): # single set of lat-lons matches station, infill from station
                    df.loc[df['elevation'] == elev_value, 'elevation_eraqc'] = 4 # see era_qaqc_flag_meanings.csv
                    df.loc[df['elevation'] == elev_value, 'elevation'] = df['elevation'].iloc[0]

                else: # lat-lon of incorrectly coded elevation does not match station lat-lon (has shifted), infill from dem
                    dem_elev_value = _grab_dem_elev_m(off_lats[0], off_lons[0])
                    df.loc[df['elevation'] == elev_value, 'elevation_eraqc'] = 3 # see era_qaqc_flag_meanings.csv
                    df.loc[df['elevation'] == elev_value, 'elevation'] = float(dem_elev_value)
            
            else: # multple pairs of lat-lon for incorrectly zero coded elevs 
                for ilat in off_lats:
                    for ilon in off_lons:
                        dem_elev_value = _grab_dem_elev_m(ilat, ilon)
                        df.loc[(df['lat'] == ilat) & (df['lon'] == ilon), 'elevation_eraqc'] = 3 # see era_qaqc_flag_meanings.csv
                        df.loc[(df['lat'] == ilat) & (df['lon'] == ilon), 'elevation'] = float(dem_elev_value)
    
    print('Elevation values post-infilling/correcting: {}'.format(df['elevation'].unique())) # testing
    print('Elevation qaqc values post-infilling/correcting: {}'.format(df['elevation_eraqc'].unique())) # testing

    # then check for in range
    # if value is present but outside of reasonable value range

    # death valley is 282 feet (85.9 m) below sea level
    # denali is ~6190 m

    if (df['elevation'].values.any() < -86.0) or (df['elevation'].values.any() > 6200.0):
        df = pd.DataFrame() # returns empty df to flag that it does not pass

    # elevation value is present and within reasonable value range
    else:
        df = df

    return df


## Time conversions
## Need function to calculate sub-hourly to hourly -- later on?

#----------------------------------------------------------------------
## Part 2 functions (individual variable/timestamp)

## NDBC and MARITIME only

def spurious_buoy_check(station, df, qc_vars):
    """
    Checks the end date on specific buoys to confirm disestablishment/drifting dates of coverage.
    If station reports data past disestablishment date, data records are flagged as suspect.
    """
    known_issues = ['NDBC_46023', 'NDBC_46045', 'NDBC_46051', 'NDBC_46044', 'MARITIME_PTAC1', 'MARITIME_PTWW1', 'MARITIME_MTYC1', 'MARITIME_MEYC1',
                    'MARITIME_SMOC1', 'MARITIME_ICAC1']
    potential_issues = ['NDBC_46290', 'NDBC_46404', 'NDBC_46212', 'NDBC_46216', 'NDBC_46220', 'NDBC_46226', 'NDBC_46227', 'NDBC_46228', 
                        'NDBC_46230', 'NDBC_46234', 'NDBC_46245', 'NDBC_46250']
                        
    if station in known_issues:
        print('{0} is flagged as suspect, checking data coverage'.format(station)) # testing

        # buoys with "data" past their disestablishment dates
        if station == 'NDBC_46023': # disestablished 9/8/2010
            for i in range(df.shape[0]):
                for j in qc_vars:
                    if (df.index[i][-1].date() >= datetime.date(2010, 9, 9)) == True:
                        df.loc[df.index[i], j] = 2 # see era_qaqc_flag_meanings.csv
            
        elif station == "NDBC_46045": # disestablished 11/1997
            for i in range(df.shape[0]):
                for j in qc_vars:
                    if (df.index[i][-1].date() >= datetime.date(1997, 12, 1)) == True:
                        df.loc[df.index[i], j] = 2 # see era_qaqc_flag_meanings.csv

        elif station == "NDBC_46051": # disestablished 4/1996, and out of range of DEM (past coastal range) but reports nan elevation
            # qaqc_elev_infill sets elevation to be assumed 0.0m, but flagging 
            for i in range(df.shape[0]):
                for j in qc_vars:
                    if (df.index[i][-1].date() >= datetime.date(1996, 5, 1)) == True:
                        df.loc[df.index[i], j] = 2 # see era_qaqc_flag_meanings.csv 
          
        elif station == "MARITIME_PTAC1": # data currently available 1984-2012, but disestablished 2/9/2022
            # only flag if new data is added after 2022 in a new data pull
            for i in range(df.shape[0]):
                for j in qc_vars:
                    if (df.index[i][-1].date() >= datetime.date(2022, 2, 9)) == True:
                        df.loc[df.index[i], j] = 2 # see era_qaqc_flag_meanings.csv

        # adrift buoy that reports valid data during adrift period (5/2/2015 1040Z to 5/3/2015 1600Z)
        elif station == "NDBC_46044":
            for i in range(df.shape[0]):
                for j in qc_vars:
                    if (df.index[i][-1] >= datetime.datetime(2015, 5, 2, 10, 40, 0)) and (df.index[i][-1] <= datetime.datetime(2015, 5, 3, 15, 50, 0)):
                        df.loc[df.index[i], j] = 2 # see era_qaqc_flag_meanings.csv

        # other known issues
        elif station == "MARITIME_PTWW1": # wind data obstructed by ferries docking at pier during day hours
            # only wind vars need flag during "day" hours, currently set for 6am to 8pm every day
            for i in range(df.shape[0]):
                if (df.index[i][-1].time() >= datetime.time(6, 0)) and (df.index[i][-1].time() <= datetime.time(20, 0)):
                    df.loc[df.index[i], "sfcWind_eraqc"] = 1 # see era_qaqc_flag_meanings.csv
                    df.loc[df.index[i], "sfcWind_dir_eraqc"] = 1 

        # elif station == "MARITIME_MTYC1" or station == "MARITIME_MEYC1": # buoy was renamed, no relocation; MTYC1 2005-2016, MEYC1 2016-2021
        #     # modify attribute/naming with note
        #     # this will get flagged in station proximity tests

        # elif station == "MARITIME_SMOC1" or station == "MARITIME_ICAC1": # buoy was renamed, small relocation (see notes); SMOC1 2005-2010, ICAC1 2010-2021
        #     # modify attribute/naming with note
        #     # this will get flagged in station proximity tests
        return df

    elif station in potential_issues: 
        # other stations have partial coverage of their full data records as well as disestablishment dates
        # if new data is added in the future, needs a manual check and added to known issue list if requires handling
        # most of these should be caught by not having a cleaned data file to begin with, so if this print statement occurs it means new raw data was cleaned and added to 2_clean_wx/
        print("{0} has a reported disestablishment date, requires manual confirmation of dates of coverage".format(station))
        for i in range(df.shape[0]):
            for j in qc_vars:
                df.loc[df.index[i], j] = 2  # see era_qaqc_flag_meanings.csv
    
        return df

    else: # station is not suspicious, move on
        return df
    

## logic check: precip accumulation amounts balance for time period
def qaqc_precip_logic_accum_amounts(df):
    """
    Ensures that precipitation accumulation amounts are consistent with reporting time frame.
    Only needs to be applied when 2 or more precipitation duration specific
    variables are present (pr_5min, pr_1h, pr_24h)
    For example: pr_5min should not be larger than pr_1h
    """
    # pr: Precipitation accumulated since last record
    # pr_5min: Precipitation accumulated in last 5 minutes
    # pr_1h: Precipitation accumulated in last hour
    # pr_24h: Precipitation accumulated from last 24 hours
    # pr_localmid: Precipitation accumulated from local midnight
        
    # rules
    # pr_5min < pr_1h < pr_24h
    # pr_localmid should never exceed pr_24h

    # determine which precipitation vars are present
    pr_vars = [col for col in df.columns if 'pr_' in col] # excludes 'pr' variable
    pr_vars = [item for item in pr_vars if "qc" not in item] # excludes raw/eraqc variable
    pr_vars = [item for item in pr_vars if "duration" not in item] # excludes duration variable (if provided)

    if len(pr_vars) == 0: # if station does not report any precipitation values, bypass
        print('station does not report a precipitation duration variable - bypassing precip logic check') # testing
        df = df

    elif len(pr_vars) == 1: # no need for amount check
        print('station does not report multiple precipitation duration variables - bypassing precip logic check') # testing
        df = df
        
    elif len(pr_vars) >= 1: 
        # checks accumulated precip vars against each other
        # noting that these flags are essentially identical in operation
        # flag assignment is logically dependent on the first var to determine which flag is placed (i.e. to determine if too larges/small)
        if 'pr_5min' in pr_vars:
            if 'pr_1h' in pr_vars:
                df.loc[df['pr_5min'] > df['pr_1h'], 'pr_5min_eraqc'] = 15 # see era_qaqc_flag_meanings.csv
            if 'pr_24h' in pr_vars:
                df.loc[df['pr_5min'] > df['pr_24h'], 'pr_5min_eraqc'] = 15 # see era_qaqc_flag_meanings.csv 
            print('Precip 5min eraqc flags (any other value than nan is an active flag!): {}'.format(df['pr_5min_eraqc'].unique())) # testing

        if 'pr_1h' in pr_vars:
            if 'pr_5min' in pr_vars:
                df.loc[df['pr_1h'] < df['pr_5min'], 'pr_1h_eraqc'] = 16 # see era_qaqc_flag_meanings.csv
            if 'pr_24h' in pr_vars:
                df.loc[df['pr_1h'] > df['pr_24h'], 'pr_1h_eraqc'] = 15 # see era_qaqc_flag_meanings.csv   
            print('Precip 1h eraqc flags (any other value than nan is an active flag!): {}'.format(df['pr_1h_eraqc'].unique())) # testing

        if 'pr_24h' in pr_vars:
            if 'pr_5min' in pr_vars:
                df.loc[df['pr_24h'] < df['pr_5min'], 'pr_24h_eraqc'] = 16 # see era_qaqc_flag_meanings.csv
            if 'pr_1h' in pr_vars:
                df.loc[df['pr_24h'] < df['pr_1h'], 'pr_24h_eraqc'] = 16 # see era_qaqc_flag_meanings.csv 
            # checks pr_24h against pr_localmid, catches an issue in pr_24h
            if 'pr_localmid' in pr_vars:
                df.loc[df['pr_24h'] < df['pr_localmid'], 'pr_24h_eraqc'] = 17 # see era_qaqc_flag_meanings.csv
            print('Precip 24h eraqc flags (any other value than nan is an active flag!): {}'.format(df['pr_24h_eraqc'].unique())) # testing

    return df


## missing value check: double check that all missing value observations are converted to NA before QA/QC
def qaqc_missing_vals(df):
    '''
    Checks data to be qaqc'ed for any errant missing values that made it through cleaning
    Converts those missing values to NAs
    Searches for missing values in 2_clean_data/missing_data_flags.csv
    '''

    missing_vals = pd.read_csv('missing_data_flags.csv')
    
    all_vars = [col for col in df.columns if 'qc' not in col]
    obs_vars = [var for var in all_vars if var not in ['lon','lat']]
    
    for item in obs_vars:
        # pull missing values which are appropriate for the range of real values for each variable 
        missing_codes = missing_vals.loc[missing_vals['variable'].str.contains(item) | missing_vals['variable'].str.contains('all')]
        
        # values in column that == missing_flag values, replace with NAs
        # note numerical vals converted to strings first to match missing_flag formatting
        df[item] = np.where(df[item].astype(str).isin(missing_codes['missing_flag']), float('NaN'), df[item])
        
        print(item)
        
    return df

## logic check: precip does not have any negative values
def qaqc_precip_logic_nonegvals(df):
    """
    Ensures that precipitation values are positive. Negative values are flagged as impossible.
    Provides handling for the multiple precipitation variables presently in the cleaned data. 
    """
    # pr_24h: Precipitation accumulated from last 24 hours
    # pr_localmid: Precipitation accumulated from local midnight
    # pr: Precipitation accumulated since last record
    # pr_1h: Precipitation accumulated in last hour
    # pr_5min: Precipitation accumulated in last 5 minutes

    # identify which precipitation vars are reported by a station
    all_pr_vars = [col for col in df.columns if 'pr' in col] # can be variable length depending if there is a raw qc var
    pr_vars = [var for var in all_pr_vars if 'qc' not in var] # remove all qc variables so they do not also run through: raw, eraqc, qaqc_process
    pr_vars = [var for var in pr_vars if 'method' not in var]
    pr_vars = [var for var in pr_vars if 'duration' not in var]

    if len(pr_vars) != 0: # precipitation variable(s) is present
        for item in pr_vars:
            print('Precip range: ', df[item].min(), '-', df[item].max()) # testing
            if (df[item] < 0).any() == True:
                df.loc[df[item] < 0, item+'_eraqc'] = 10 # see era_qaqc_flag_meanings.csv

            print('Precipitation eraqc flags (any other value than nan is an active flag!): {}'.format(df[item+'_eraqc'].unique())) # testing

    else: # station does not report precipitation
        print('station does not report precipitation - bypassing precip logic check') # testing
        df = df

    return df

  
## sensor height - air temperature
def qaqc_sensor_height_t(xr_ds, file_to_qaqc):
    '''
    Checks if temperature sensor height is within 2 meters above surface +/- 1/3 meter tolerance.
    If missing or outside range, temperature value for station is flagged to not proceed through QA/QC.
    '''
    
    # Check if thermometer height is missing
    if (np.isnan(xr_ds.thermometer_height_m)):
        file_to_qaqc['tas_eraqc'] = file_to_qaqc['tas_eraqc'].fillna(6) # see era_qaqc_flag_meanings.csv
    
    else: # sensor height present
        # Check if thermometer height is within 2 m +/- 1/3 m
        if(xr_ds.thermometer_height_m >= (2 - 1/3) and xr_ds.thermometer_height_m <= (2 + 1/3)):
            file_to_qaqc = file_to_qaqc
                
        else: 
            # Thermometer height present but outside 2m +/- tolerance
            file_to_qaqc['tas_eraqc'] = file_to_qaqc['tas_eraqc'].fillna(7)
            
    return file_to_qaqc

## sensor height - wind
def qaqc_sensor_height_w(xr_ds, file_to_qaqc):
    '''
    Checks if wind sensor height is within 10 meters above surface +/- 1/3 meter tolerance.
    If missing or outside range, wind speed and direction values for station are flagged to not proceed through QA/QC.
    '''
        
    # Check if anemometer height is missing
    if np.isnan(xr_ds.anemometer_height_m):
        file_to_qaqc['sfcWind_eraqc'] = file_to_qaqc['sfcWind_eraqc'].fillna(8) # see era_qaqc_flag_meanings.csv
        file_to_qaqc['sfcWind_dir_eraqc'] = file_to_qaqc['sfcWind_dir_eraqc'].fillna(8)
    
    else: # sensor height present
        if xr_ds.anemometer_height_m >= (10 - 1/3) and xr_ds.anemometer_height_m <= (10 + 1/3):
            # Check if anemometer height is within 10 m +/- 1/3 m
            file_to_qaqc = file_to_qaqc
                    
        else: 
            # Anemometer height present but outside 10m +/- tolerance
            file_to_qaqc['sfcWind_eraqc'] = file_to_qaqc['sfcWind_eraqc'].fillna(9)
            file_to_qaqc['sfcWind_dir_eraqc'] = file_to_qaqc['sfcWind_dir_eraqc'].fillna(9)
                
    return file_to_qaqc

## flag values outside world records for North America
# temp, dewpoint, windspeed, sea level pressure
def qaqc_world_record(df):
    '''
    Checks if temperature, dewpoint, windspeed, or sea level pressure are outside North American world records
    If outside minimum or maximum records, flags values
    '''
    
    # world records from HadISD protocol, cross-checked with WMO database
    # https://wmo.asu.edu/content/world-meteorological-organization-global-weather-climate-extremes-archive
    T_X = {"North_America":329.92} #K
    T_N = {"North_America":210.15} #K
    D_X = {"North_America":329.85} #K
    D_N = {"North_America":173.15} #K
    W_X = {"North_America":113.2} #m/s
    W_N = {"North_America":0.} #m/s
    S_X = {"North_America":108330} #Pa
    S_N = {"North_America":87000} #Pa

    maxes = {"tas": T_X, "tdps": D_X, "tdps_derived": D_X, "sfcWind": W_X, "psl": S_X}
    mins = {"tas": T_N, "tdps": D_N, "tdps_derived": D_N, "sfcWind": W_N, "psl": S_N}
    
    # column names to check against world record limits
    wr_cols = ['tas', 'tdps_derived', 'tdps', 'sfcWind', 'psl']

    # subset data to variables to check
    wr_vars = [col for col in df.columns if col in wr_cols]

    for item in wr_vars:
        if ((df[item] < mins[item]['North_America']) | (df[item] > maxes[item]['North_America'])).any() == True:
                df.loc[(df[item] < mins[item]['North_America']) | (df[item] > maxes[item]['North_America']), item+'_eraqc'] = 11
    
    return df


## cross-variable logic checks
# dew point must not exceed air temperature
def qaqc_crossvar_logic_tdps_to_tas(df):
    """
    Checks that dewpoint temperature does not exceed air temperature.
    If fails, only dewpoint temperature is flagged.
    """ 

    # First check that tdps and/or tdps_derived are provided
    dew_vars = [col for col in df.columns if 'tdps' in col]
    all_dew_vars = [var for var in dew_vars if 'qc' not in var] # remove all qc variables so they do not also run through: raw, eraqc

    if len(all_dew_vars) != 0: # dew point is present
        for dew_var in all_dew_vars:            
            # check values for physical constraint
            df.loc[df[dew_var] > df['tas'], dew_var+"_eraqc"] = 12  # see qaqc_flag_meanings.csv

            print('{0} eraqc flags (any other value than nan is an active flag!): {1}'.format(dew_var, df[dew_var+'_eraqc'].unique())) # testing

    else: # station does not report dew point temperature
        print('station does not report dew point temperature - bypassing temperature cross-variable logic check') # testing
        df = df

    return df

# wind direction must be 0 if wind speed is 0
def qaqc_crossvar_logic_calm_wind_dir(df):
    """
    Checks that wind direction is zero when wind speed is also zero.
    If fails, wind direction is flagged. # only flag wind direction?
    """

    # Noting that a wind direction value of 0 is a valid value
    # Only a problem when wind speed is also 0, where 0 now means no winds for there to be a direction

    # First check that wind direction is provided
    if 'sfcWind_dir' in df.columns:
        # First, identify calm winds but with incorrect wind directions
        df.loc[(df['sfcWind'] == 0) & # calm winds
            (df['sfcWind_dir'] != 0) & # direction is not 0
            (df['sfcWind_dir'].isnull() == False), # exclude directions that are null/nan
                'sfcWind_dir_eraqc'] = 13 # see qaqc_flag_meanings.csv

        # Next, identify non-zero winds but with incorrect wind directions
        # Non-zero northerly winds should be coded as 360 deg, not 0 deg
        df.loc[(df['sfcWind'] != 0) & # non-calm winds
            (df['sfcWind_dir'] == 0) & # direction is zero
            (df['sfcWind_dir'].isnull() == False), # exclude directions that are null/nan
            'sfcWind_dir_eraqc'] = 14 # see era_qaqc_flag_meanings.csv

        print('sfcWind_dir eraqc flags (any value other than nan is an active flag!): {}'.format(df['sfcWind_dir_eraqc'].unique()))

    else: # station does not report wind direction
        print('station does not report wind direction - bypassing wind cross-variable logic check') # testing
        df = df

    return df



## distribution checks
# flag unusually frequent values - if any one value has more than 50% of data in the bin
def create_bins(data, bin_size=0.25):
    '''Create bins from data covering entire data range'''

    # set up bins
    b_min = np.floor(np.nanmin(data))
    b_max = np.ceil(np.nanmax(data))
    bins = np.arange(b_min - bin_size, b_max + (3. * bin_size), bin_size)

    return bins

def synergistic_flag(df, num_temp_vars):  
    '''
    In frequent values, if air temp is flagged, dew point is also flagged, and vice versa.
    Applies appropriate flag in corresponding vars
    '''
    if 'tas' in num_temp_vars and 'tdps' in num_temp_vars:
        df.loc[df['tas_eraqc'] == 22, 'tdps_eraqc'] = 22
        df.loc[df['tdps_eraqc'] == 22, 'tas_eraqc'] = 22

    if 'tas' in num_temp_vars and 'tdps_derived' in num_temp_vars:
        df.loc[df['tas_eraqc'] == 22, 'tdps_derived_eraqc'] = 22
        df.loc[df['tdps_derived_eraqc'] == 22, 'tas_eraqc'] = 22    

    if 'tas' in num_temp_vars and 'tdps' in num_temp_vars and 'tdps_derived' in num_temp_vars:
        df.loc[df['tas_eraqc'] == 22, 'tdps_eraqc'] = 22
        df.loc[df['tdps_eraqc'] == 22, 'tas_eraqc'] = 22
        df.loc[df['tas_eraqc'] == 22, 'tdps_derived_eraqc'] = 22
        df.loc[df['tdps_derived_eraqc'] == 22, 'tas_eraqc'] = 22
        df.loc[df['tdps_eraqc'] == 22, 'tdps_derived_eraqc'] = 22
        df.loc[df['tdps_derived_eraqc'] == 22, 'tdps_eraqc'] = 22
            
    return df

def qaqc_frequent_vals(df, plots=True):
    '''Frequent values check:
        - Initially > 50% of all data in current 1 degC/hPa bin 
        - out of "this and +/- 3 bins for all data to highlight with >30 (obs?) in the bin
        - On yearly basis using highlighted bins with 50% of data and >=20 obs in this and +/- 3 bins OR
        - 90% data and >=10 observations in this and +/-3 bins
        - for seasons, bin size thresholds are reduced to 20, 15, and 10 respectively
        
        Note: tas and tdps are synergistic
            - if t is bad, tdps is also removed, and vice versa
    '''
    
    # this check is only done on air temp, dewpoint temp, and pressure
    vars_to_remove = ['qc', 'duration', 'method']
    vars_to_include = ['tas', 'tdps', 'ps'] # list of var substrings to remove if present in var
    vars_to_check = [var for var in df.columns if any(True for item in vars_to_include if item in var) and not any(True for item in vars_to_remove if item in var)]

    df = df.reset_index() 
    df['month'] = pd.to_datetime(df['time']).dt.month # sets month to new variable
    df['year'] = pd.to_datetime(df['time']).dt.year # sets year to new variable
    
    for var in vars_to_check:
        # set-up flag vars
        df[var+'_eraqc'] = np.nan
        
        # first scans suspect values using entire record
        # all years
        df = frequent_bincheck(df, var, data_group='all')

        # if no values are flagged as suspect, end function, no need to proceed
        if len(df.loc[df[var+'_eraqc'] == 100]) == 0:
            print('No unusually frequent values detected for entire {} observation record'.format(var))
            # goes to seasonal check, no bypass

        else:
            # year by year
            # then scans for each value on a year-by-year basis to flag if they are a problem within that year
                # DECISION: the annual check uses the unfiltered data
                # previously flagged values are included here -- this would interfere with our entire workflow
            df = frequent_bincheck(df, var, data_group='annual')

        # seasonal scan (JF+D, MAM, JJA, SON) 
        # each season is scanned over entire record to identify problem values
        # only flags applied on annual basis using the three months on their own
        # NOTE: HadISD approach is to use the current year's december, rather than the preceeding december

        # seasonal version because seasonal shift in distribution of temps/dewpoints can reveal hidden values
        # all years
        df = frequent_bincheck(df, var, data_group='seasonal_all') ## DECISION: December is from the current year
        if len(df.loc[df[var+'_eraqc'] == 100]) == 0:
            print('No unusually frequent values detected for seasonal {} observation record'.format(var))
            continue # bypasses to next variable

        else:
            print('Unusually frequent values detected in seasonal distribution, continuining to annual check')
            # year by year --> December selection must be specific
            df = frequent_bincheck(df, var, data_group='seasonal_annual')    
                
                
        # remove any lingering preliminary flags, data passed check
        df.loc[df[var+'_eraqc'] == 100, var+'_eraqc'] = np.nan
        
    # synergistic flag on tas and tdps/tdps_derived
    # first establish at least tas and one tdps var present
    temp_vars = ['tas', 'tdps', 'tdps_derived']
    num_temp_vars = [var for var in vars_to_check if var in temp_vars]
    if len(num_temp_vars) != 1 and 'tas' in num_temp_vars:
        # proceed to synergistic check
        df = synergistic_flag(df, num_temp_vars)
    
    # plots item
    if plots==True:
        for var in vars_to_check:
            if 22 in df[var+'_eraqc'].values: # only plot a figure if a value is flagged
                # histogram
                frequent_vals_plot(df, var)

                # entire timeseries figure
                flagged_timeseries_plot(df, flag_to_viz=frequent_flags)
        
    return df


def frequent_vals_plot(df, var):
    '''
    Produces a histogram of the diagnostic histogram per variable, 
    and any bin that is indicated as "too frequent" by the qaqc_frequent_vals test 
    is visually flagged
    ''' 
    bins = create_bins(df[var], 1)
    ax = df.plot.hist(column=var, bins=bins, alpha=0.5)
    
    # plot flagged values
    
    # first identify which values are flagged
    vals_to_flag = df.loc[df[var+'_eraqc'] == 22][var].unique()
    bars_to_flag = []
    for i in vals_to_flag:
        if math.isnan(i) == False:
            bars_to_flag.append(math.floor(i))
            
    # flag bars if too frequent
    for bar in ax.patches:
        x = bar.get_x() + 0.5 * bar.get_width()
        if x+0.5 in bars_to_flag: # right tail
            bar.set_color('r')

    # plot aesthetics
    plt.title('Frequent value check: {}'.format(df['station'].unique()[0]),
             fontsize=10);
    
    
    # save figure to AWS
    network = df['station'].unique()[0].split('_')[0]
    
    bucket_name = 'wecc-historical-wx'
    directory = '3_qaqc_wx'
    img_data = BytesIO()
    plt.savefig(img_data, format='png')
    img_data.seek(0)

    s3 = boto3.resource('s3')
    bucket = s3.Bucket(bucket_name)
    figname = 'qaqc_frequent_value_check_{0}_{1}'.format(df['station'].unique()[0], var)
    bucket.put_object(Body=img_data, ContentType='image/png',
                     Key='{0}/{1}/qaqc_figs/{2}.png'.format(
                     directory, network, figname))

    # close figures to save memory
    plt.close()

def frequent_bincheck(df, var, data_group):
    '''Approach: 
        - histograms created with 0.5 or 1.0 or hpa increments (depending on accuracy of instrument)
        - each bin compared to the three on either side
        - if this bin contains more than half the total population of the seven bins combined
        - and more than 30 observations over the station record (20 for seasonal)
        - then histogram bin is highlighted for further investigation
        - minimum number limit imposted to avoid removing true tails of distribution
    '''    
    
    # seasons
    szns = [[3,4,5], [6,7,8], [9,10,11], [12,1,2]] 
    
    # bin sizes: using 1 degC for tas/tdps, and 1 hPa for ps vars
    ps_vars = ['ps', 'ps_altimeter', 'psl']
    
    ## TEMPORARY BUG FIX ON PSL UNIT ==================================================================
    if len(str(df.loc[df.index == df['psl'].first_valid_index(), 'psl'].values[0]).split('.')[0]) <= 4:
        df['psl'] = df['psl'] * 100
    ## END TEMPORARY FIX ==============================================================================
    
    if var in ps_vars: 
        bin_s = 100 # all of our pressure vars are in Pa, convert to 100 Pa bin size
    else:
        bin_s = 1 
    
    # all data/annual checks
    if data_group == 'all':
        bins = create_bins(df[var], bin_size=bin_s) 
        bar_counts, bins = np.histogram(df[var], bins=bins)
        flagged_bins = bins_to_flag(bins, bar_counts)
        
        # flag values in that bin as suspect
        if len(flagged_bins) != 0:
            for sus_bin in flagged_bins:
                # indicate as suspect bins
                    # DECISION: preliminary flag? and then remove if okay/reset to nan?
                df.loc[(df[var]>=sus_bin) & (df[var]<=sus_bin+1), 
                       var+'_eraqc'] = 100 # highlight for further review flag, either overwritten with real flag or removed in next step
    
    #============================================================================================================
       
    elif data_group == 'annual':
        for yr in df.year.unique():
            df_yr = df.loc[df['year'] == yr]
            bins = create_bins(df_yr[var], bin_size=bin_s) # using 1 degC/hPa bin width
            bar_counts, bins = np.histogram(df_yr[var], bins=bins)
            flagged_bins = bins_to_flag(df_yr, bar_counts, bin_main_thresh=20, secondary_bin_main_thresh=10)
            
            if len(flagged_bins) != 0:
                for sus_bin in flagged_bins:
                    print('Flagging bin: ', sus_bin)
                    df.loc[(df['year']==yr) & (df[var]>=sus_bin) & (df[var]<=sus_bin+1), 
                           var+'_eraqc'] = 22 # see era_qaqc_flag_meanings.csv
    
    
    #============================================================================================================
    # seasonal checks require special handling
    elif data_group == 'seasonal_all':
        for szn in szns:
            df_szn = df.loc[(df['month']==szn[0]) | (df['month']==szn[1]) | (df['month']==szn[2])]
            bins = create_bins(df_szn[var], bin_size=bin_s) # using 1 degC/hPa bin width
            bar_counts, bins = np.histogram(df_szn[var], bins=bins)
            flagged_bins = bins_to_flag(df_szn[var], bar_counts, bin_main_thresh=20, secondary_bin_main_thresh=20)
            
            if len(flagged_bins) != 0:
                for sus_bin in flagged_bins:
                    df.loc[((df['month']==szn[0]) | (df['month']==szn[1]) | (df['month']==szn[2])) & 
                           (df[var]>=sus_bin) & (df[var]<=sus_bin+1),
                           var+'_eraqc'] = 100 # highlight for further review flag, either overwritten with real flag or removed in next step
                    
                    
    #============================================================================================================
                
    elif data_group == 'seasonal_annual':        
        for yr in df.year.unique():
            for szn in szns:
                # all seasons except winter
                if szn != [12,1,2]:
                    df_szn = df.loc[(df['year']==yr) & 
                                    ((df['month']==szn[0]) | (df['month']==szn[1]) | (df['month']==szn[2]))]                    
                    
                    if yr==df.loc[df.index[-1],'year']:
                        if len(df_szn)==0:
                            break # after last season in last year
                    
                    bins = create_bins(df_szn[var], bin_size=bin_s) # using 1 degC/hPa bin width
                    bar_counts, bins = np.histogram(df_szn[var], bins=bins)
                    flagged_bins = bins_to_flag(df_szn[var], bar_counts, bin_main_thresh=15, secondary_bin_main_thresh=10)
                    
                    if len(flagged_bins) != 0:
                        for sus_bin in flagged_bins:
                            print('Flagging bin: ', sus_bin)
                            df.loc[(df['year']==yr) & 
                                  ((df['month']==szn[0]) | (df['month']==szn[1]) | (df['month']==szn[2])) &
                                   (df[var]>=sus_bin) & (df[var]<=sus_bin+1),
                                  var+'_eraqc'] = 22 # see era_qaqc_flag_meanings.csv

                # special handling for winter because of december
                else:
                    df_yr = df.loc[df['year'] == yr] # that year's jan, feb, and wrong dec            
                    df_jf = df_yr.loc[df['month'] != 12] # that specific year's jan and feb

                    df_d = df.loc[(df['year'] == yr-1) & (df['month'] == 12)] # previous year's dec
                    if len(df_d) == 0: # catching very first year instance
                        df_djf = df_jf 
                        print('Winter season: proceeding with just Jan/Feb, no previous Dec') ## DECISION

                    else:
                        print('Winter season: concatenating previous Dec')
                        df_djf = pd.concat([df_d, df_jf])
                    
                    bins = create_bins(df_djf[var], bin_size=bin_s) # using 1 degC/hPa bin width
                    bar_counts, bins = np.histogram(df_djf[var], bins=bins)
                    flagged_bins = bins_to_flag(df_djf[var], bar_counts, bin_main_thresh=15, secondary_bin_main_thresh=10)

                    if len(flagged_bins) != 0:
                        for sus_bin in flagged_bins:
                            print('Flagging bin: ', sus_bin)
                            # flag jan feb
                            df.loc[(df['year']==yr) & 
                                   ((df['month']==szn[1]) | (df['month']==szn[2])) &
                                   ((df[var]>=sus_bin) & (df[var]<=sus_bin+1)),
                                  var+'_eraqc'] = 22 # see era_qaqc_flag_meanings.csv
                            # flag correct dec
                            df.loc[((df['year']==yr-1) & (df['month']==szn[0])) &
                                   ((df[var]>=sus_bin) & (df[var]<=sus_bin+1)),
                                   var+'_eraqc'] = 22 # see era_qaqc_flag_meanings.csv
                
    return df


def bins_to_flag(bins, bar_counts, bin_main_thresh=30, secondary_bin_main_thresh=30):
    '''Returns the specific bins to flag as suspect'''
    bins_to_flag = [] # list of bins that will be flagged
    
    for i in range(0, len(bar_counts)):
        # identify main bin + 3 on either side
        bin_start = i-3
        bin_end = i+4

        # need handling for first 3 blocks as there is no front
        if i < 3:
            bin_start = 0

        bin_block_sum = bar_counts[bin_start:bin_end].sum() # num of obs in the 7-bin block
        bin_main_sum = bar_counts[i] # num of obs in main bin

        # determine whether main bin is more than half sum in 7-block bin
        bin_block_50 = bin_block_sum * 0.5 # primary check at 50%
        bin_block_90 = bin_block_sum * 0.9 # secondary check at 90%

        if (bin_main_sum > bin_block_50) == True: 
            # ensure that bin_main_sum is greater than bin_main_thresh
            if bin_main_sum > bin_main_thresh:
                bins_to_flag.append(math.floor(bins.values[i]))
                
                # annual/seasonal check
                if (bin_main_sum > bin_block_90) == True:
                    if bin_main_sum > secondary_bin_main_thresh:
                        bins_to_flag.append(math.floor(bins.values[i])) 
                
            else: # less than bin_main_thresh obs in bin_main_sum, do not indicate as suspect
                continue
                
    return bins_to_flag # returns a list of values that are suspicious   

#----------------------------------------------------------------------
# To do
# establish false positive rate
# unit tests on each of these functions(?)
