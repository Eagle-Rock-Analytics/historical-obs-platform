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

# flag unusual gaps within the monthly distribution bins
def qaqc_dist_gaps(df, plots=False):
    '''
    Identifies if there are any unusual gaps in the monthly distribution for any variable.
    Flags if there is a gap in the distribution, and outputs and saves a figure when data is flagged. 
    '''

    # run through every var, excluding qaqc/duration/method vars
    vars_to_remove = ['qc', 'duration', 'method', 'lat', 'lon', 'elevation'] # list of var substrings to exclude if present in var
    vars_to_check = [var for var in df.columns if not any(True for item in vars_to_remove if item in var)] # remove all non-primary variables

    # first need to set a minimum threshold for # of obs to build distribution check

    print('Checking for gaps in monthly distribution for: {}'.format(var))

    # if data is flagged, print statement and save figure
    print('Unusual gap in monthly distribution identified for {0} in month of {1} - flagged and figure saved for analysis'.format(var, MONTH))

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
    
#----------------------------------------------------------------------
# To do
# establish false positive rate
# unit tests on each of these functions(?)
