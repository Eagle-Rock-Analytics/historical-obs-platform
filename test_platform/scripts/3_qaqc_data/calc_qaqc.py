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
import shapely
import xarray as xr
from scipy import stats
import matplotlib.pyplot as plt
import boto3
from io import BytesIO

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


#======================================================================
## PART 1 functions (whole station/network)

#----------------------------------------------------------------------
# missing spatial coords (lat-lon)
def qaqc_missing_latlon(df, verbose=True):
    """
    Checks if latitude and longitude is missing for a station.
    If missing, station is flagged to not proceed through QA/QC.
    """

    # latitude or longitude
    variables = list(df.columns)
    if "lon" not in variables or "lat" not in variables:
        return None

    if df['lat'].isnull().all():
        return None

    if df['lon'].isnull().all():
        return None

    # df['lon'] = df['lon'].fillna(method="pad")
    # df['lat'] = df['lon'].fillna(method="pad")
    
    return df
        
#----------------------------------------------------------------------
# in bounds of WECC
def qaqc_within_wecc(df, verbose=True):
    """
    Checks if station is within terrestrial & marine WECC boundaries.
    If outside of boundaries, station is flagged to not proceed through QA/QC.
    """

    t = gp.read_file(wecc_terr).iloc[0].geometry  ## Read in terrestrial WECC shapefile.
    m = gp.read_file(wecc_mar).iloc[0].geometry   ## Read in marine WECC shapefile.
    pxy = shapely.geometry.Point(df['lon'].mean(), df['lat'].mean())
    if pxy.within(t) or pxy.within(m):
        return df
    else:
        return None

#----------------------------------------------------------------------
# elevation
def _grab_dem_elev_m(lats_to_check, lons_to_check, verbose=True):
    """
    Pulls elevation value from the USGS Elevation Point Query Service, 
    lat lon must be in decimal degrees (which it is after cleaning)
    Modified from: 
    https://gis.stackexchange.com/questions/338392/getting-elevation-for-multiple-lat-long-coordinates-in-python
    """
    url = r'https://epqs.nationalmap.gov/v1/json?'

    dem_elev_short = np.ones_like(lats_to_check)*np.nan
    
    for i,lat,lon in zip(range(len(lats_to_check)), lats_to_check, lons_to_check):
        # define rest query params
        params = {
            'output': 'json',
            'x': lon,
            'y': lat,
            'units': 'Meters'
        }

        # format query string and return value
        result = requests.get((url + urllib.parse.urlencode(params)))
        dem_elev_long = float(result.json()['value'])
        # make sure to round off lat-lon values so they are not improbably precise for our needs
        # dem_elev_short[i] = '{:.2f}'.format(dem_elev_long) 
        dem_elev_short[i] = np.round(dem_elev_long, decimals=2) 

    return dem_elev_short.astype("float")

#----------------------------------------------------------------------
def qaqc_elev_infill(df, verbose=True):
    """
    Checks if elevation is NA/missing. If missing, fill in elevation from either DEM or station.
    Some stations have all nan elevation values (e.g., NDBC, MARITIME)
    Some stations have single/few but not all nan elevation values (e.g., otherisd, asosawos)
    """    
    if verbose:
        print('Elevation values pre-infilling: {}'.format(df['elevation'].unique()))
        print('Elevation eraqc values pre-infilling: {}'.format(df['elevation_eraqc'].unique())) # testing

    # elev values missing
    isNan = df['elevation'].isnull()

    # if all are missing
    if isNan.any():
        if isNan.all():
            dem_elev_values = _grab_dem_elev_m([df['lat'].iloc[0]], 
                                               [df['lon'].iloc[0]],
                                               verbose=verbose)    
            df.loc[:, 'elevation'] = dem_elev_values[0]    
            df.loc[:, 'elevation_eraqc'] = 3    
        else:
            # if some missing
            try:
                if df['lon'].is_unique and df['lat'].is_unique:
                    dem_elev_values = _grab_dem_elev_m([df['lat'].iloc[0]], 
                                                       [df['lon'].iloc[0]],
                                                       verbose=verbose)
                    df.loc[:, 'elevation'] = dem_elev_values[0]
                    df.loc[:, 'elevation_eraqc'] = 3
                else:
                    dem_elev_values = _grab_dem_elev_m([df.loc[isNan, 'lat']], 
                                                       [df.loc[isNan, 'lon']],
                                                        verbose=verbose)
                    df.loc[isNan, 'elevation'] = dem_elev_values
                    df.loc[isNan, 'elevation_eraqc'] = 3
                return df
            
            # elevation cannot be obtained from DEM
            except:
                if verbose:
                    print("Elevation cannot be obtained from DEM")
                return None
    else:
        return df

#----------------------------------------------------------------------
def qaqc_elev_range(df, verbose=True):
    """
    Checks valid values to identify suspicious elevation values that are larger than 10m in difference
    Checks if valid elevation value is outside of range of reasonable values for WECC region.
    If outside range, station is flagged to not proceed through QA/QC.
    """
    # first check for suspicious values
    # elev_vals = df['elevation'].unique() # identify how many different elevation "values" are present

    # elevation values flagged as incorrectly coded
    # uses a threshold of 10m different from the station elevation to identify suspicious elevations
    # control = _grab_dem_elev_m([df.loc['lat'].iloc[0]], 
    #                            [df.loc['lon'].iloc[0]])[0]

    # TO DO:
    # This is the original version, but what about if the first value is bad? (look up to the DEM version)
    control = df['elevation'].iloc[0]
    isOff = np.logical_or(df['elevation'] > control + 10,
                          df['elevation'] < control - 10)
    if isOff.any():
        # in-fill if value is missing
        try:
            if df['lon'].is_unique and df['lat'].is_unique:
                dem_elev_values = _grab_dem_elev_m([df['lat'].iloc[0]], 
                                                   [df['lon'].iloc[0]])
                df.loc[ifOff, 'elevation'] = dem_elev_values[0]
                df.loc[isOff, 'elevation_eraqc'] = 3
            else:
                dem_elev_values = _grab_dem_elev_m([df.loc[ifOff, 'lat']], 
                                                   [df.loc[isOff, 'lon']])
                df.loc[ifOff, 'elevation'] = dem_elev_values
                df.loc[isOff, 'elevation_eraqc'] = 3

        # elevation cannot be obtained from DEM
        except:
            return None
    else:
        return df

    if verbose:
        print('Elevation values post-infilling/correcting: {}'.format(df['elevation'].unique())) # testing
        print('Elevation qaqc values post-infilling/correcting: {}'.format(df['elevation_eraqc'].unique())) # testing
    
    # then check for in range if value is present but outside of reasonable value range
    # death valley is 282 feet (85.9 m) below sea level, denali is ~6190 m
    
    isOut = df['elevation'].iloc[0] < -86.0 or df['elevation'].iloc[0]>6200.0
    if isOut.any():
        return None
    else:
        return df
    
## Time conversions
## Need function to calculate sub-hourly to hourly -- later on?

#======================================================================
## Part 2: Variable logic checks

#-----------------------------------------------------------------------------
## logic check: precip does not have any negative values
def qaqc_precip_logic_nonegvals(df, verbose=True):
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
    all_pr_vars = [var for var in df.columns if 'pr' in var] # can be variable length depending if there is a raw qc var
    pr_vars = [var for var in all_pr_vars if 'qc' not in var] # remove all qc variables so they do not also run through: raw, eraqc, qaqc_process
    pr_vars = [var for var in pr_vars if 'method' not in var]
    pr_vars = [var for var in pr_vars if 'duration' not in var]

    if not pr_vars: # precipitation variable(s) is not present
        print('station does not report precipitation - bypassing precip logic check')
        return None
    else:
        for item in pr_vars:
            if verbose:
                print('Precip range: ', df[item].min(), '-', df[item].max()) # testing
            isNeg = df[item] < 0
            df.loc[isNeg, item+'_eraqc'] = 10 # see era_qaqc_flag_meanings.csv

            if verbose:
                print('Precipitation eraqc flags (any other value than nan is an active flag!):' + 
                      '{}'.format(df[item+'_eraqc'].unique())) # testing
    return df

#----------------------------------------------------------------------
## logic check: precip accumulation amounts balance for time period
def qaqc_precip_logic_accum_amounts(df, verbose=True):
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

    # identify which precipitation vars are reported by a station
    all_pr_vars = [var for var in df.columns if 'pr' in var] # can be variable length depending if there is a raw qc var
    pr_vars = [var for var in all_pr_vars if 'qc' not in var] # remove all qc variables so they do not also run through: raw, eraqc, qaqc_process
    pr_vars = [var for var in pr_vars if 'method' not in var]
    pr_vars = [var for var in pr_vars if 'duration' not in var]
    
    # if station does not report any precipitation values, or only one, bypass
    if len(pr_vars) == 0 or len(pr_vars) == 1:
        return df

    # checks accumulated precip vars against each other
    # noting that these flags are essentially identical in operation
    # flag assignment is logically dependent on the first var to determine 
    # which flag is placed (i.e. to determine if too larges/small)

    # checks accumulated precip vars against each other
    # noting that these flags are essentially identical in operation
    # flag assignment is logically dependent on the first var to determine which flag is placed (i.e. to determine if too larges/small)
    
    #:::::
    if 'pr_5min' in pr_vars:
        if 'pr_1h' in pr_vars:
            isBad = df['pr_5min'] > df['pr_1h']
            df.loc[isBad, 'pr_5min_eraqc'] = 15 # see era_qaqc_flag_meanings.csv
            
        if 'pr_24h' in pr_vars:
            isBad = df['pr_5min'] > df['pr_24h']
            df.loc[isBad, 'pr_5min_eraqc'] = 15 # see era_qaqc_flag_meanings.csv
        if verbose:
            print('Precip 5min eraqc flags (any other value than nan is an active flag!):' + 
                  '{}'.format(df['pr_5min_eraqc'].unique())) # testing

    #:::::
    if 'pr_1h' in pr_vars:
        if 'pr_5min' in pr_vars:
            isBad = df['pr_1h'] < df['pr_5min']
            df.loc[isBad, 'pr_1h_eraqc'] = 16 # see era_qaqc_flag_meanings.csv
            
        if 'pr_24h' in pr_vars:
            isBad = df['pr_1h'] > df['pr_24h']
            df.loc[isBad, 'pr_1h_eraqc'] = 15 # see era_qaqc_flag_meanings.csv
        if verbose:
            print('Precip 1h eraqc flags (any other value than nan is an active flag!):' + 
                  '{}'.format(df['pr_1h_eraqc'].unique())) # testing

    #:::::
    if 'pr_24h' in pr_vars:
        if 'pr_5min' in pr_vars:
            isBad = df['pr_24h'] < df['pr_5min']
            df.loc[isBad, 'pr_24h_eraqc'] = 16 # see era_qaqc_flag_meanings.csv
        
        if 'pr_1h' in pr_vars:
            isBad = df['pr_24h'] < df['pr_1h']
            df.loc[isBad, 'pr_24h_eraqc'] = 16 # see era_qaqc_flag_meanings.csv        
        
        if 'pr_localmid' in pr_vars:
            isBad = df['pr_24h'] < df['pr_localmid']
            df.loc[isBad, 'pr_24h_eraqc'] = 17 # see era_qaqc_flag_meanings.csv
            
        if verbose:
            print('Precip 24h eraqc flags (any other value than nan is an active flag!):' + 
                  '{}'.format(np.unique(df['pr_24h_eraqc']))) # testing

    return df

#======================================================================
## Part 3 functions (individual variable/timestamp)
## NDBC and MARITIME only

#----------------------------------------------------------------------
def spurious_buoy_check(df, qc_vars, verbose=True):
    """
    Checks the end date on specific buoys to confirm disestablishment/drifting dates of coverage.
    If station reports data past disestablishment date, data records are flagged as suspect.
    """
    known_issues = ['NDBC_46023', 'NDBC_46045', 'NDBC_46051', 'NDBC_46044', 'MARITIME_PTAC1', 'MARITIME_PTWW1', 'MARITIME_MTYC1', 'MARITIME_MEYC1',
                    'MARITIME_SMOC1', 'MARITIME_ICAC1']
    potential_issues = ['NDBC_46290', 'NDBC_46404', 'NDBC_46212', 'NDBC_46216', 'NDBC_46220', 'NDBC_46226', 'NDBC_46227', 'NDBC_46228', 
                        'NDBC_46230', 'NDBC_46234', 'NDBC_46245', 'NDBC_46250']

    # remove elevation_qc var from remainder of analyses so it does not also get flagged -- 
    # confirm with final qaqc process
    if "elevation_eraqc" in qc_vars:
        qc_vars.remove("elevation_eraqc") 
    
    # Extract station name from dataset encoding
    station = df.index[0][0]
    
    if station in known_issues:
        if verbose:
            print('{0} is flagged as suspect, checking data coverage'.format(station)) # testing
        
        # buoys with "data" past their disestablishment dates
        if station == 'NDBC_46023': # disestablished 9/8/2010
            isBad = df['time'] >= np.datetime64("2010-09-09")
            for new_var in qc_vars:
                # # Retrieve originaldd var name
                # var = new_var.split("_eraqc")[0]
                df.loc[isBad, new_var] = 2 # see era_qaqc_flag_meanings.csv
            
        elif station == "NDBC_46045": # disestablished 11/1997
            isBad = df['time'] >= np.datetime64("1997-12-01")
            for new_var in qc_vars:
                df.loc[isBad, new_var] = 2 # see era_qaqc_flag_meanings.csv

        elif station == "NDBC_46051": # disestablished 4/1996, and out of range of DEM (past coastal range) but reports nan elevation
            isBad = df['time'] >= np.datetime64("1996-05-01")
            for new_var in qc_vars:
                df.loc[isBad, new_var] = 2 # see era_qaqc_flag_meanings.csv

        elif station == "MARITIME_PTAC1": # data currently available 1984-2012, but disestablished 2/9/2022
            # only flag if new data is added after 2022 in a new data pull
            isBad = df['time'] >= np.datetime64("2022-02-09")
            for new_var in qc_vars:
                df.loc[isBad, new_var] = 2 # see era_qaqc_flag_meanings.csv

        # adrift buoy that reports valid data during adrift period (5/2/2015 1040Z to 5/3/2015 1600Z)
        elif station == "NDBC_46044":
            isBad = df['time'] >= np.datetime64("2015-05-02 10:40:00") and df['time'] <= np.datetime64("2015-05-03 15:50:00")
            for new_var in qc_vars:
                df.loc[isBad, new_var] = 2 # see era_qaqc_flag_meanings.csv
                
        # other known issues
        elif station == "MARITIME_PTWW1": # wind data obstructed by ferries docking at pier during day hours
            # only wind vars need flag during "day" hours, currently set for 6am to 8pm every day
            isBad = df['time'] >= np.datetime64("1900-01-01 06:00:00") and df['time'] <= np.datetime64("1900-01-01 20:00:00")
            
            df.loc[isBad, "sfcWind_eraqc"] = 1
            df.loc[isBad, "sfcWind_dir_eraqc"] = 1 # see era_qaqc_flag_meanings.csv

        # elif station == "MARITIME_MTYC1" or station == "MARITIME_MEYC1": # buoy was renamed, no relocation; MTYC1 2005-2016, MEYC1 2016-2021
        #     # modify attribute/naming with note
        #     # this will get flagged in station proximity tests

        # elif station == "MARITIME_SMOC1" or station == "MARITIME_ICAC1": # buoy was renamed, small relocation (see notes); SMOC1 2005-2010, ICAC1 2010-2021
        #     # modify attribute/naming with note
        #     # this will get flagged in station proximity tests

    elif station in potential_issues: 
        # other stations have partial coverage of their full data records as well as disestablishment dates
        # if new data is added in the future, needs a manual check and added to known issue list if requires handling
        # most of these should be caught by not having a cleaned data file to begin with, so if this print statement occurs it means new raw data was cleaned and added to 2_clean_wx/
        if verbose:
            print("{0} has a reported disestablishment date, requires manual confirmation of dates of coverage".format(station))
        
        for new_var in qc_vars:
            df.loc[:, new_var] = 2 # see era_qaqc_flag_meanings.csv

    return df

#-----------------------------------------------------------------------------
## sensor height - wind
def qaqc_sensor_height_w(df, verbose=True):
    '''
    Checks if wind sensor height is within 10 meters above surface +/- 1/3 meter tolerance.
    If missing or outside range, wind speed and direction values for station are flagged to not proceed through QA/QC.
    '''
    # try:
    if True:
        # Check if anemometer height is missing
        isHeightMissing = df['anemometer_height_m'].isnull().any()

        if isHeightMissing:
            # df.loc[:,'sfcWind_eraqc'] = 8 # see era_qaqc_flag_meanings.csv
            # df.loc[:,'sfcWind_dir_eraqc'] = 8
            df['sfcWind_eraqc'] = 8 # see era_qaqc_flag_meanings.csv
            df['sfcWind_dir_eraqc'] = 8

        else: # sensor height present
            # Check if anemometer height is within 10 m +/- 1/3 m
            isHeightWithin = df['anemometer_height_m'][0] >= (10 - 1/3) and df['anemometer_height_m'][0] <= (10 + 1/3)
            # Anemometer height present but outside 10m +/- tolerance
            if not isHeightWithin:
                df['sfcWind_eraqc'] = 9
                df['sfcWind_dir_eraqc'] = 9 
        return df

    else:
    # except Exception as e:
        if verbose:
            print("qaqc_sensor_height_w failed with Exception: {}".format(e))
        return None

#-----------------------------------------------------------------------------
## sensor height - air temperature
def qaqc_sensor_height_t(df, verbose=True):
    '''
    Checks if temperature sensor height is within 2 meters above surface +/- 1/3 meter tolerance.
    If missing or outside range, temperature value for station is flagged to not proceed through QA/QC.
    '''
    try:
        # Check if thermometer height is missing   
        isHeightMissing = df['thermometer_height_m'].isnull().any()

        if isHeightMissing:
            # df.loc[:,'tas_eraqc'] = 6 # see era_qaqc_flag_meanings.csv
            df['tas_eraqc'] = 6 # see era_qaqc_flag_meanings.csv
        else:
            isHeightWithin = np.logical_and(df['thermometer_height_m'] >= (2 - 1/3),
                                            df['thermometer_height_m'] <= (2 + 1/3))

            # Thermometer height present but outside 10m +/- tolerance
            if not isHeightWithin:
                # df.loc[:, 'tas_eraqc'] = 7
                df['tas_eraqc'] = 7

        return df
    except Exception as e:
        if verbose:
            print("qaqc_sensor_height_w failed with Exception: {}".format(e))
        return None

#-----------------------------------------------------------------------------------------
## flag values outside world records for North America
# temp, dewpoint, windspeed, sea level pressure
def qaqc_world_record(df, verbose=True):
    '''
    Checks if temperature, dewpoint, windspeed, or sea level pressure are outside North American world records
    If outside minimum or maximum records, flags values
    '''
    try:
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

        # variable names to check against world record limits
        wr_vars = ['tas', 'tdps_derived', 'tdps', 'sfcWind', 'psl']

        for var in wr_vars:
            if var in list(df.columns):
                isOffRecord = np.logical_or(df[var] < mins[var]['North_America'],
                                            df[var] > maxes[var]['North_America'])
                if isOffRecord.any():
                    df.loc[isOffRecord, var + '_eraqc'] = 11
        return df
    except Exception as e:
        if verbose:
            print("qaqc_world_record failed with Exception: {}".format(e))
        return None

#-----------------------------------------------------------------------------------------
## cross-variable logic checks
# dew point must not exceed air temperature
def qaqc_crossvar_logic_tdps_to_tas(df, verbose=True):
    """
    Checks that dewpoint temperature does not exceed air temperature.
    If fails, only dewpoint temperature is flagged.
    """ 
    try:
        # First check that tdps and/or tdps_derived are provided
        dew_vars = [col for col in df.columns if 'tdps' in col]
        all_dew_vars = [var for var in dew_vars if 'qc' not in var] # remove all qc variables so they do not also run through: raw, eraqc

        # dew point is not present
        if not all_dew_vars:
            if verbose:
                print('station does not report dew point temperature - bypassing temperature cross-variable logic check')
        # dew point is present
        else:
            for var in all_dew_vars: 
                isBad = df[var] > df['tas']
                df.loc[isBad, var + '_eraqc'] = 12 # see qaqc_flag_meanings.csv
                if verbose:
                    print('{0} eraqc flags (any other value than nan is an active flag!): {1}'.
                          format(var, df[var + '_eraqc'].unique()))
        return df
    
    except Exception as e:
        print("qaqc_crossvar_logic_tdps_to_tas failed with Exception: {}".format(e))
        return None

#-----------------------------------------------------------------------------------------
# wind direction must be 0 if wind speed is 0
def qaqc_crossvar_logic_calm_wind_dir(df, verbose=True):
    """
    Checks that wind direction is zero when wind speed is also zero.
    If fails, wind direction is flagged. # only flag wind direction?
    """
    try:
        # Noting that a wind direction value of 0 is a valid value
        # Only a problem when wind speed is also 0, where 0 now means no winds for there to be a direction

        # First check that wind direction is provided
        if 'sfcWind_dir' not in df.columns:
            if verbose:
                print('station does not report wind direction - bypassing wind cross-variable logic check')
                return df
            
        # First, identify calm winds but with incorrect wind directions
        isCalm = df['sfcWind'] == 0
        isDirNotZero = df['sfcWind_dir'] != 0
        isNotNan = ~df['sfcWind_dir'].isnull()
        isBad = isCalm & isDirNotZero & isNotNan
        
        df.loc[isBad, 'sfcWind_dir_eraqc'] = 13 # see qaqc_flag_meanings.csv
        
        # Next, identify non-zero winds but with incorrect wind directions
        # Non-zero northerly winds should be coded as 360 deg, not 0 deg
        isNotCalm = df['sfcWind'] != 0
        isDirZero = df['sfcWind_dir'] == 0
        isBad = isNotCalm & isDirZero & isNotNan
        
        df.loc[isBad, 'sfcWind_dir_eraqc'] = 14 # see qaqc_flag_meanings.csv
        
        if verbose:
            print('sfcWind_dir eraqc flags (any value other than nan is an active flag!): {}'.
                  format(df['sfcWind_dir_eraqc'].unique()))
        return df
        
    except Exception as e:
        print("qaqc_crossvar_logic_calm_wind_dir failed with Exception: {}".format(e))
        return None
    
#-----------------------------------------------------------------------------------------
# flag unusual gaps within the monthly distribution bins
def qaqc_dist_gaps_part1(df):
    """
    Part 1 / monthly check
        - compare anomalies of monthly median values
        - standardize against interquartile range
        - compare stepwise from the middle of the distribution outwards
        - asymmetries are identified and flagged if severe
    Goal: identifies individual suspect observations and flags the entire month  
    """
    
    # in order to grab the time information more easily -- would prefer not to do this
    df = df.reset_index() 
    df['month'] = pd.to_datetime(df['time']).dt.month # sets month to new variable
    df['year'] = pd.to_datetime(df['time']).dt.year # sets year to new variable
    
    # calculate monthly medians
    df_anom = df.sub(df.resample('M', on='time').transform('median', numeric_only=True))
    df_anom['time'] = df['time'] # add time column back in to do quantiles

    # standardize against calendar-monIQR range
    df_q1 = df_anom.resample('M', on='time').transform('quantile', 0.25, numeric_only=True)
    df_q3 = df_anom.resample('M', on='time').transform('quantile', 0.75, numeric_only=True)
    df_iqr = df_q3 - df_q1
    df_anom_iqr = df_anom / df_iqr
    
    # run through every var, excluding qaqc/duration/method vars
    vars_to_remove = ['index','station','qc', 'duration', 'method', 'lat', 'lon', 'elevation', 'time', 'month', 'year'] # list of var substrings to exclude if present in var
    vars_to_check = [var for var in df.columns if not any(True for item in vars_to_remove if item in var)] # remove all non-primary variables
    
    for var in vars_to_check:
        # add _eraqc column for each variable
        df[var+'_eraqc'] = np.nan # default value of nan
        
        # "inflated to 4°C or hPa for those months with very small IQR"
        # accounts for any seasonal cycle in variance
        small_iqr_var_check = ['tas', 'tdps', 'tdps_derived', 'ps', 'psl', 'psl_altimeter', 'ps_derived']
        if var in small_iqr_var_check:
            if (np.abs(df_anom_iqr[var].max()) + np.abs(df_anom_iqr[var].min())) < 4:
                print('small var check') # testing for occurrence 
                df_anom_iqr[var] = np.linspace(-2, 2, len(df)) # unsure this is the correct way to do this - come back

        # standardized anomalies are ranked (necessary?) and calculate median
        std_med = df_anom_iqr.median() # will be 0 if inflated to range of 4

        # add standardized anomaly median to IQR-standardized data
        df_std_med = df_anom_iqr + std_med
        df_std_med['time'] = df['time'] # add time columns back in... again
        df_std_med['year'] = df['year']
        df_std_med['month'] = df['month']

        # identify where any obs are +/- 5 IQR away from standardized anomaly median
        if len(df_std_med.loc[np.abs(df_std_med[var]) > 5]) != 0:

            bad_idxs = df_std_med.loc[np.abs(df_std_med[var]) > 5].index.tolist() # grab indices of suspect obs
            print('{} suspicious {} observations present, flagging appropriate months'.format(len(bad_idxs), var))

            for i in bad_idxs:
                bad_yr = df.iloc[df.index == i]['year'].values[0]
                bad_mon = df.iloc[df.index == i]['month'].values[0]
                print('Flagging: {0}/{1}'.format(bad_mon, bad_yr))

                # identify all indices for months encapsulating suspect obs
                bad_obs_per_month = df.loc[(df['year'] == bad_yr) & (df['month'] == bad_mon)]
                all_idx_to_flag = bad_obs_per_month.index

                for i in all_idx_to_flag: # flag all indices in those months
                    df.loc[df.index == i, var+'_eraqc'] = 18 # see era_qaqc_flag_meanings.csv # DOUBLE CHECK VALUE

        else:
            print('Part 1: PASS. All {} obs are within +/- 5 IQR range'.format(var))

    return df

#===================================================================================================
## Part 4: Distribution functions
## distribution gap function helpers

#-----------------------------------------------------------------------------------------

def monthly_med(df):
    """Part 1: Calculates the monthly median"""
    return df.resample('M', on='time').median(numeric_only=True)

#-----------------------------------------------------------------------------------------

def iqr_range(df, month, var):
    """Part 1: Calculates the monthly interquartile range"""
    q1 = df.groupby('month').quantile(0.25, numeric_only=True)
    q3 = df.groupby('month').quantile(0.75, numeric_only=True)
    iqr_df = q3 - q1
    
    iqr_val = iqr_df.loc[iqr_df.index == month]
    
    # # inflated to 4°C or 4 hPa for months with very small IQR
    # var_check = ['tas', 'tdps', 'tdps_derived', 'ps', 'psl', 'psl_altimeter']
    # if iqr_val[var].values < 4:
    #     if var in var_check:
    #         iqr_val[var].values = 4
    
    return iqr_val[var].values

#-----------------------------------------------------------------------------------------

def standardize_iqr(df, var):
    """Part 2: Standardizes data against the interquartile range

    Returns:
        array
    """
    q1 = df[var].quantile(0.25)
    q3 = df[var].quantile(0.75)
    iqr = q3 - q1

    return (df[var].values - df[var].median()) / iqr
    
#-----------------------------------------------------------------------------------------

def median_clim(df, month, var):
    '''Part 2: Calculate climatological median for a specific month and variable'''

    clim = df[var].median(numeric_only=True)

    return clim

#-----------------------------------------------------------------------------------------

def standardized_anom(df, month, var):
    """
    Part 1: Calculates the monthly anomalies standardized by IQR range
    
    Returns:
        arr_std_anom: array of monthly standardized anomalies for var
    """
    
    df_monthly_med = monthly_med(df)
    df_clim_med = clim_med(df)
    
    arr_anom = (df_monthly_med.loc[df_monthly_med['month'] == month][var].values -
                df_clim_med.loc[df_clim_med.index == month][var].values)
        
    arr_std_anom = arr_anom / iqr_range(df, month, var)
    
    return arr_std_anom
    
#-----------------------------------------------------------------------------------------

def standardized_median_bounds(df, month, var, iqr_thresh=5):
    """Part 1: Calculates the standardized median"""
    std_med = df.loc[df['month'] == month][var].median() # climatological median for that month
    
    lower_bnd = std_med - (iqr_thresh * iqr_range(df, month, var))
    upper_bnd = std_med + (iqr_thresh * iqr_range(df, month, var))
    
    return (std_med, lower_bnd[0], upper_bnd[0])
    
#-----------------------------------------------------------------------------------------

def qaqc_dist_whole_stn_bypass_check(df, vars_to_check, min_num_months=5):
    """Part 1: Checks the number of valid observation months in order to proceed through monthly distribution checks. Identifies whether a station record has too 
    few months and produces a fail pass flag. 
    """
    
    # set up df
    df = df.reset_index() 
    df['month'] = pd.to_datetime(df['time']).dt.month # sets month to new variable
    df['year'] = pd.to_datetime(df['time']).dt.year # sets year to new variable
             
    # set up a "pass_flag" to determine if station proceeds through distribution function
    pass_flag = 'pass'
    
    for var in vars_to_check:
        # add _eraqc column for each variable
        df[var+'_eraqc'] = np.nan # default value of nan

        for month in range(1,13):

            # first check num of months in order to continue
            month_to_check = df.loc[df['month'] == month]

            # check for number of obs years
            if (len(month_to_check.year.unique()) < 5):
                df[var+'_eraqc'] = 18 # see era_qaqc_flag_meanings.csv
                pass_flag = 'fail'

    err_statement = '{} has too short of an observation record to proceed through the monthly distribution qa/qc checks -- bypassing station'.format(
                    df['station'].unique()[0])
    
    if pass_flag == 'fail':
        print(err_statement)
                
    return (df, pass_flag) 

#-----------------------------------------------------------------------------------------

def qaqc_dist_var_bypass_check(df, vars_to_check, min_num_months=5):
    """
    Part 1: Checks the number of valid observation months per variable to proceed through monthly distribution checks.
    Primarily assesses whether if null values persist for a month
    """
        
    for var in vars_to_check:
        for month in range(1,13):
            monthly_df = df.loc[df['month']==month]
            
            # if all values are null for that month across years
            if monthly_df[var].isnull().all() == True:
                df[var+'_eraqc'] = 19 # see era_qaqc_flag_meanings.csv
            
            # if not all months have nans, need to assess how many years do
            elif monthly_med(df).loc[monthly_med(df)['month'] == month][var].isna().sum() > min_num_months:
                print('{} has less than {} months of valid {} observations -- bypassing check'.format(
                df['station'].unique()[0], min_num_months, var))
                
                df[var+'_eraqc'] = 19 # see era_qaqc_flag_meanings.csv
        
    return df

#-----------------------------------------------------------------------------------------

## distribution gap plotting helpers
def create_bins(data, bin_size=0.25):
    '''Create bins from data covering entire data range'''

    # set up bins
    b_min = np.floor(np.nanmin(data))
    b_max = np.ceil(np.nanmax(data))
    bins = np.arange(b_min - bin_size, b_max + (3. * bin_size), bin_size)

    return bins

#-----------------------------------------------------------------------------------------

def _plot_format_helper(var):
    """Helper function for plots"""

    pr_vars = ['pr', 'pr_5min', 'pr_1h', 'pr_24h', 'pr_localmid']
    ps_vars = ['ps', 'psl', 'psl_altimeter']
    
    if var == 'tas':
        ylab = 'Air Temperature at 2m'
        unit = 'K'
        
    elif var == 'tdps' or var == 'tdps_derived':
        ylab = 'Dewpoint Temperature'
        unit = 'K'
        
    elif var == 'sfcWind':
        ylab = 'Surface Wind Speed'
        unit = '${m s^-1}$'
        
    elif var == 'sfcWind_dir':
        ylab = 'Surface Wind Direction'
        unit = 'degrees'
        
    elif var == 'rsds':
        ylab = 'Surface Radiation'
        unit = '${W m^-2}$'

    elif var == 'hurs':
        ylab = 'Humidity'
        unit = '%'
        
    elif var in pr_vars:
        ylab = 'Precipitation' # should be which precip var it is
        unit = 'mm'

    elif var in ps_vars:
        ylab = 'Pressure' # should eventually be what pressure var it is
        unit = 'Pa'
        
    return (ylab, unit)

#-----------------------------------------------------------------------------------------

## distribution gap flagging functions
def qaqc_dist_gap_part1(df, iqr_thresh=5, plot=True):
    """
    Part 1 / monthly check
        - compare anomalies of monthly median values
        - standardize against interquartile range
        - compare stepwise from the middle of the distribution outwards
        - asymmetries are identified and flagged if severe
    Goal: identifies suspect months and flags all obs within month
    
    Note: PRELIMINARY: This function has not been fully evaluated or finalized in full qaqc process. Thresholds/decisions may change with refinement.
        - iqr_thresh preliminarily set to 5 years, pending revision
    """
    
    # run through every var, excluding qaqc/duration/method vars
    vars_to_remove = ['index','station','qc','duration','method','lat','lon','elevation','time','month','year','sfcWind_dir','hurs'] # list of var substrings to exclude if present in var
    vars_to_check = [var for var in df.columns if not any(True for item in vars_to_remove if item in var)] # remove all non-primary variables
    
    for var in vars_to_check:
        for month in range(1,13): 

            # per variable bypass check
            df = qaqc_dist_var_bypass_check(df, vars_to_check) # flag here is 19
            if 19 in df[var+'_eraqc']:
                continue # skip variable 

            # station has above min_num_months number of valid observations, proceed with dist gap check
            else:
                # calculate monthly climatological median, and bounds
                mid, low, high = standardized_median_bounds(df, month, var, iqr_thresh=iqr_thresh)

                # calculate monthly median per month
                df_month = monthly_med(df)

                for i in df_month.loc[df_month['month'] == month][var]:
                    if (i < low) or (i > high):
                        year_to_flag = (df_month.loc[(df_month[var]==i) & 
                                           (df_month['month']==month)]['year'].values[0])
                        print('Median {} value for {}-{} is beyond the {}*IQR limits -- flagging month'.format(
                            var,
                            month, 
                            int(year_to_flag),
                            iqr_thresh)
                        )

                        # flag all obs in that month
                        df.loc[(df['month']==month) & 
                               (df['year']==year_to_flag), var+'_eraqc'] = 20 # see era_qaqc_flag_meanings.csv

        if plot==True:
            for month in range(1,13):
                for var in vars_to_check:
                    if 19 not in df[var+'_eraqc'].values: # don't plot a figure if it's all nans/not enough months
                        dist_gap_part1_plot(df, month, var, flagval=20, iqr_thresh=iqr_thresh,
                                            network=df['station'].unique()[0].split('_')[0])
                
    return df

#-----------------------------------------------------------------------------------------

def qaqc_dist_gap_part2(df, plot=True):
    """
    Part 2 / monthly check
        - compare all obs in a single month, all years
        - histogram created from all obs and gaussian distribution fitted
        - threshold values determined using positions where fitted freq falls below y=0.1
        - rounds outwards to next integer plus one
        - going outwards from center, distribution is scanned for gaps which occur outside threshold
        - obs beyond gap are flagged
    Goal: identifies individual suspect observations and flags the entire month  
    """
    # run through every var, excluding qaqc/duration/method vars
    vars_to_remove = ['index','station','qc','duration','method','lat','lon','elevation','time','month','year','sfcWind_dir','hurs'] # list of var substrings to exclude if present in var
    vars_to_check = [var for var in df.columns if not any(True for item in vars_to_remove if item in var)] # remove all non-primary variables
        
    # whole station bypass check first
    df, pass_flag = qaqc_dist_whole_stn_bypass_check(df, vars_to_check)
    
    if pass_flag != 'fail':
        
        for var in vars_to_check:
            for month in range(1,13):
                
                # per variable bypass check
                df = qaqc_dist_var_bypass_check(df, vars_to_check) # flag here is 19
                if 19 in df[var+'_eraqc']:
                    continue # skip variable 
                
                # station has above min_num_months number of valid observations, proceed with dist gap check
                else:
                    # from center of distribution, scan for gaps (where bin = 0)
                    # when gap is found, and it is at least 2x bin width
                    # any bins beyond end of gap + beyond threshold value are flagged
                    
                    # subset by month
                    df = df.loc[df['month'] == month]
                    
                    # standardize against IQR range
                    df_month_iqr = iqr_standardize(df, var)

                    # determine number of bins
                    bins = create_bins(df_month_iqr)
                    
                    # pdf
                    mu = np.nanmean(df_month_iqr)
                    sigma = np.nanstd(df_month_iqr)

                    y, left_bnd, right_bnd = pdf_bounds(df_month_iqr, mu, sigma, bins)
                    
                    # identify gaps as below y=0.1 from histogram, not pdf                    
                    y_hist, bins = np.histogram(df_iqr, bins=bins, density=True)
                    
                    # identify climatology and iqr baselines in order to flag
                    iqr_baseline = iqr_range(df, month=month, var=var)
                    clim = median_clim(df, month=month, var=var)
                                        
                    # gaps are only flagged for values beyond left_bnd, right_bnd, as long as gap is 2*bin_width (2*0.25)
                    # considering that the # of bins for threshold is (4,7) from y=0.1
                    # safe to assume that gap is present if values >0.1 outside of left_bnd, right_bnd
                    bins_beyond_left_bnd = np.argwhere(bins <= left_bnd)
                    if len(bins_beyond_left_bnd) != 0: 
                        for data in bins_beyond_left_bnd:
                            if y_hist[data] > 0.1: # bins with data > 0.1 beyond left_bnd
                                
                                # identify values beyond left bnd
                                vals_to_flag = clim + (left_bnd * iqr_baseline) # left_bnd is negative
                                df.loc[df[var] <= vals_to_flag[0], var+'_eraqc'] = 21 # see era_qaqc_flag_meanings.csv


                    bins_beyond_right_bnd = np.argwhere(bins >= right_bnd)
                    if len(bins_beyond_right_bnd) != 0:
                        for data in bins_beyond_right_bnd:
                            if y_hist[data] > 0.1: # bins with data > 0.1 beyond right_bnd
                                
                                # identify values beyond right bnd
                                vals_to_flag = clim + (right_bnd * iqr_baseline) # upper limit threshold
                                df.loc[df[var] >= vals_to_flag[0], var+'_eraqc'] = 21 # see era_qaqc_flag_meanings.csv
                    
    ## Question: Do we need "all", "flagged_only", "none" options instead?
    if plot==True:
        for month in range(1,13):
            for var in vars_to_check:
                if 19 not in df[var+'_eraqc'].values: # don't plot a figure if it's all nans/not enough months
                    dist_gap_part2_plot(df, month, var,
                                        network=df['station'].unique()[0].split('_')[0])
    
    return df  

#-----------------------------------------------------------------------------------------

## distribution gap plotting functions
def dist_gap_part1_plot(df, month, var, flagval, iqr_thresh, network):
    '''Produces a timeseries plots of specific months and variables for part 1 of the unusual gaps function.
    Any variable that is flagged is noted'''

    # grab data by months
    df = df.loc[df['month'] == month]
        
    # grab flagged data
    flag_vals = df.loc[df[var + '_eraqc'] == flagval]
    
    # plot valid data
    ax = df.plot.scatter(x='time', y=var, label='Pass')
    
    # plot flagged data
    flag_vals.plot.scatter(ax=ax, x='time', y=var, color='r', label='Flagged')
    # should be consistent with other plots - I like Hector's open circles around flagged values

    # plot climatological median and threshold * IQR range
    mid, low_bnd, high_bnd = standardized_median_bounds(df, month, var, iqr_thresh=5)
    
    plt.axhline(y=mid, color='k', lw=0.5, label='Climatological monthly median')
    plt.fill_between(x=df['time'],
                    y1=low_bnd,
                    y2=high_bnd,
                    alpha=0.25, color='0.75', 
                    label='{} * IQR range'.format(iqr_thresh))
    
    # plot aesthetics
    plt.legend(loc='best')
    ylab = _plot_format_helper(var)
    plt.ylabel('{} [{}]'.format(ylab[0], ylab[1]));
    plt.xlabel('')
    plt.title('Distribution gap check pt 1: {0} / month: {1}'.format(
        df['station'].unique()[0],
        month), 
              fontsize=10);
    
    # save to AWS
    bucket_name = 'wecc-historical-wx'
    directory = '3_qaqc_wx'
    img_data = BytesIO()
    plt.savefig(img_data, format='png')
    img_data.seek(0)
    
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(bucket_name)
    figname = 'qaqc_dist_gap_check_part1_{0}_{1}_{2}'.format(df['station'].unique()[0], var, month)
    bucket.put_object(Body=img_data, ContentType='image/png',
                 Key='{0}/{1}/qaqc_figs/{2}.png'.format(
                 directory, network, figname))
    
    # close figures to save memory
    plt.close()

#-----------------------------------------------------------------------------------------

def dist_gap_part2_plot(df, month, var, network):
    '''
    Produces a histogram of the monthly standardized distribution
    with PDF overlay and threshold lines where pdf falls below y=0.1.
    Any bin that is outside of the threshold is visually flagged
    ''' 

    # select month
    df = df.loc[df['month'] == month]
    
    # standardize against IQR range
    df_month_iqr = standardize_iqr(df, var)
    
    # determine number of bins
    bins = create_bins(df_month_iqr)
    
    # plot histogram
    ax = plt.hist(df_month_iqr, bins=bins, log=False, density=True, alpha=0.3);
    xmin, xmax = plt.xlim()
    plt.ylim(ymin=0.1)

    # plot pdf
    mu = np.nanmean(df_month_iqr)
    sigma = np.nanstd(df_month_iqr)
    y = stats.norm.pdf(bins, mu, sigma)
    l = plt.plot(bins, y, 'k--', linewidth=1)
    
    # add vertical lines to indicate thresholds where pdf y=0.1
    pdf_bounds = np.argwhere(y > 0.1)

    # find first index
    left_bnd = round(bins[pdf_bounds[0][0] -1])
    right_bnd = round(bins[pdf_bounds[-1][0] + 1])
    thresholds = (left_bnd - 1, right_bnd + 1)

    plt.axvline(thresholds[1], color='r') # right tail
    plt.axvline(thresholds[0], color='r') # left tail
    
    # flag (visually) obs that are beyond threshold
    for bar in ax[2].patches:
        x = bar.get_x() + 0.5 * bar.get_width()
        if x > thresholds[1]: # right tail
            bar.set_color('r')
        elif x < thresholds[0]: # left tail
            bar.set_color('r')

    # title and useful annotations
    plt.title('Distribution gap check, {0}: {1}'.format(df['station'].unique()[0], var), fontsize=10);
    plt.annotate('Month: {}'.format(month), xy=(0.025, 0.95), xycoords='axes fraction', fontsize=8);
    plt.annotate('Mean: {}'.format(round(mu,3)), xy=(0.025, 0.9), xycoords='axes fraction', fontsize=8);
    plt.annotate('Std.Dev: {}'.format(round(sigma,3)), xy=(0.025, 0.85), xycoords='axes fraction', fontsize=8);
    plt.ylabel('Frequency (obs)')
    
    # save figure to AWS
    bucket_name = 'wecc-historical-wx'
    directory = '3_qaqc_wx'
    img_data = BytesIO()
    plt.savefig(img_data, format='png')
    img_data.seek(0)
    
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(bucket_name)
    figname = 'qaqc_dist_gap_check_part2_{0}_{1}_{2}'.format(df['station'].unique()[0], var, month)
    bucket.put_object(Body=img_data, ContentType='image/png',
                     Key='{0}/{1}/qaqc_figs/{2}.png'.format(
                     directory, network, figname))

    # close figures to save memory
    plt.close()

#-----------------------------------------------------------------------------------------

def flagged_timeseries_plot(df, flag_to_viz):
    
    vars_to_remove = ['index','station','qc','duration','method','lat','lon','elevation','time','month','year','sfcWind_dir'] # list of var substrings to exclude if present in var
    vars_to_check = [var for var in df.columns if not any(True for item in vars_to_remove if item in var)] # remove all non-primary variables
    
    # can pass a list of flags
    for flag in flag_to_viz:
    
        # assess where each variable has flagged values
        for var in vars_to_check:
            flagged_data = df.loc[df[var+'_eraqc'] == flag]

            # only produce a plot if there is flagged values
            if len(flagged_data) == 0:
                continue

            # plot
            ax = df.plot.scatter(x='time', y=var, color='k', s=0.8, label='Valid')

            # plot flagged data
            flagged_data.plot.scatter(ax=ax, x='time', y=var, color='r', s=0.9, label='Flag: {}'.format(flag))

            # plot aesthetics
            plt.legend(loc='best', ncol=2)
            ylab = _plot_format_helper(var)
            plt.ylabel('{} [{}]'.format(ylab[0], ylab[1]));
            plt.xlabel('')
            plt.title('{0}'.format(df['station'].unique()[0]), fontsize=10);

        # save to AWS
        bucket_name = 'wecc-historical-wx'
        directory = '3_qaqc_wx'
        img_data = BytesIO()
        plt.savefig(img_data, format='png')
        img_data.seek(0)

        s3 = boto3.resource('s3')
        bucket = s3.Bucket(bucket_name)
        figname = 'flagged_timeseries_{0}_{1}'.format(df['station'].unique()[0], var)
    bucket.put_object(Body=img_data, ContentType='image/png',
                 Key='{0}/{1}/qaqc_figs/{2}.png'.format(
                 directory, network, figname))
    
    # close figures to save memory
    plt.close()

#-----------------------------------------------------------------------------------------

def qaqc_unusual_gaps(df, iqr_thresh=5, plots=True):
    '''
    Runs all parts of the unusual gaps function, with a whole station bypass check first.
    
    Note: PRELIMINARY: This function has not been fully evaluated or finalized in full qaqc process. Thresholds/decisions may change with refinement.
        - iqr_thresh preliminarily set to 5 years, pending revision
    '''

    # bypass check
    vars_to_remove = ['index','station','qc','duration','method','lat','lon','elevation','time','month','year','sfcWind_dir','hurs'] # list of var substrings to exclude if present in var
    vars_to_check = [var for var in df.columns if not any(True for item in vars_to_remove if item in var)] # remove all non-primary variables
        
    # whole station bypass check first
    df, pass_flag = qaqc_dist_whole_stn_bypass_check(df, vars_to_check)
    
    if pass_flag == 'fail':
        return df
    else:
        df_part1 = qaqc_dist_gap_part1(df, iqr_thresh, plots)
        df_part2 = qaqc_dist_gap_part2(df_part1, plots)

        if plots == True:
            for var in vars_to_check:
                if (19 not in df[var+'_eraqc'].values) and (20 in df[var+'_eraqc'].values or 21 in df[var+'_eraqc']): # don't plot a figure if it's all nans/not enough months
                    flagged_timeseries_plot(df_part2, flag_to_viz = [19, 20, 21])
    
    return df_part2

## Part 4: Distribution functions
## distribution gap function helpers
#===================================================================================================

#===================================================================================================
# UNUSUAL LARGE JUMPS DEV
#---------------------------------------------------------------------------------------------------
# Plot helper for unusual large jumps plots
def _plot_format_helper_spikes(var):
    """
    Helper function for unusual large jumps plots
    
    Input:
    -----
        var [str] : variable name being plotted
    Ouput:
    ----- ylab, unit, mins[var]['North_America'], maxes[var]['North_America'])
        ylab [str] : variable name for y-axis
        unit [str] : units of variable
        miny [float] : min var value for y axis
        maxy [float] : max var value for y axis      
    """

    pr_vars = ['pr', 'pr_5min', 'pr_1h', 'pr_24h', 'pr_localmid']
    ps_vars = ['ps', 'psl', 'psl_altimeter']
    
    if var == 'tas':
        ylab = 'Air Temperature at 2m'
        unit = 'K'
        
    elif var == 'tdps' or var == 'tdps_derived':
        ylab = 'Dewpoint Temperature'
        unit = 'K'
        
    elif var == 'sfcWind':
        ylab = 'Surface Wind Speed'
        unit = '${m s^-1}$'
        
    elif var == 'sfcWind_dir':
        ylab = 'Surface Wind Direction'
        unit = 'degrees'
        
    elif var == 'rsds':
        ylab = 'Surface Radiation'
        unit = '${W m^-2}$'
        
    elif var == 'hurs':
        ylab = 'Humidity'
        unit = '%'
        
    elif var in pr_vars:
        ylab = 'Precipitation' # should be which precip var it is
        unit = 'mm'

    elif var in ps_vars:
        ylab = 'Pressure' # should eventually be what pressure var it is
        unit = 'Pa'
        
    T_X = {"North_America":329.92} #K
    T_N = {"North_America":210.15} #K
    D_X = {"North_America":329.85} #K
    D_N = {"North_America":173.15} #K
    W_X = {"North_America":113.2}  #m/s
    W_N = {"North_America":0.}     #m/s
    S_X = {"North_America":108330} #Pa
    S_N = {"North_America":87000}  #Pa

    maxes = {"tas": T_X, "tdps": D_X, "tdps_derived": D_X, "sfcWind": W_X, "psl": S_X, "ps": S_X}
    mins =  {"tas": T_N, "tdps": D_N, "tdps_derived": D_N, "sfcWind": W_N, "psl": S_N, "ps": S_N}
    miny = mins[var]['North_America']
    maxy = maxes[var]['North_America']
    
    return ylab, unit, miny, maxy

#---------------------------------------------------------------------------------------------------
def unusual_jumps_plot(df, var, flagval=22, dpi=None, local=True):
    """
    Plots unusual large jumps qaqc result and uploads it to AWS (if local, also writes to local folder)
    Input:
    -----
        df [pd.Dataframe] : station pd.DataFrame from qaqc pipeline
        var [str] : variable name
        flagval [int] : flag value to plot (22 for unusual large jumps)
        dpi [int] : resolution for png plots
        local [bool] : if True, saves plot locally, else: only saves plot to AWS
    Ouput:
    ----- 
        None
    """
    
    # grab flagged data
    flag_vals = df.loc[df[var + '_eraqc'] == flagval]   
    
    # Create figure
    fig,ax = plt.subplots(figsize=(10,3))

    # Plot variable and flagged data
    df[var].plot(ax=ax, marker=".", ms=4, lw=1, color="k", alpha=0.5, label="Original data")
    
    flag_label = "{:.4f}% of data flagged".format(100*len(df.loc[df[var+"_eraqc"]==22, var])/len(df))
    df.loc[df[var+"_eraqc"]==22, var].plot(ax=ax, marker="o", ms=7, lw=0, mfc="none", color="C3", label=flag_label)    
    
    legend = ax.legend(loc=0, prop={'size': 8})    
        
    station = df['station'].unique()[0]
    network = station.split('_')[0]
    
    # Plot aesthetics
    ylab, units, miny, maxy = _plot_format_helper_spikes(var)
    ylab = '{} [{}]'.format(ylab, units)
    
    ax.set_ylabel(ylab)
    ax.set_xlabel('')
    
    # We can set ylim since this function is supposed to be run after other QAQC functions (including world records)
    miny = max(miny, df[var].min())
    maxy = min(maxy, df[var].max())
    ax.set_ylim(miny,maxy)
    
    title = 'Unusual large jumps check: {0}'.format(station)
    ax.set_title(title, fontsize=10)
    
    # save to AWS
    bucket_name = 'wecc-historical-wx'
    directory = '3_qaqc_wx'
    figname = 'qaqc_figs/qaqc_unusual_large_jumps_{0}_{1}'.format(station, var)
    key = '{0}/{1}/{2}.png'.format(directory, network, figname)
    img_data = BytesIO()
    fig.savefig(img_data, format='png', dpi=dpi, bbox_inches="tight")
    fig.savefig(figname+".png", format='png', dpi=dpi, bbox_inches="tight")
    img_data.seek(0)
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(bucket_name)
    bucket.put_object(Body=img_data, ContentType='image/png', Key=key)
    plt.close()
    
    return 

#---------------------------------------------------------------------------------------------------
def potential_spike_check(potential_spike, diff, crit, hours_diff):
    """
    Checks for neccessary conditions for a potential spike to be an actual spike:
     - Spikes are considered 1-value spike up to 3-values spike
     - Difference right before the spike should be lower than half the critical value
     - Difference at the actual spike must be higher than the critical value
     - Differences within the multi-value spike must lower than half the critical value
     - Difference right after (spike exit) the spike should be higher than the critical value and
       of opposite sign of the actual spike

    Input:
    -----
            potential_spike [pandas series] : bool pd.Series with True on potential spike location
            diff [pandas series] : float pd.Series with differences in the test variable
            crit [pandas series] : float pd.Series with the critical value for the differences in the test variable
            crit [pandas series] : float pd.Series with the hour differences between data points in the test variable
    Output:
    ------
            df [pandas dataframe] : input df with added `var`_spike column True where data matches the spike conditions
    """

    potential_spike = potential_spike.copy(deep=True)
    
    ind = np.where(potential_spike)[0]
    spikes = pd.Series(np.zeros_like(potential_spike).astype("bool"), index=potential_spike.index)
    dates = pd.Series(potential_spike.index.values)
    # print(len(ind))
    # print(dates[ind])
    
    for i in ind:
        
        #Ignore edges for now
        if i==1 or i>=len(potential_spike)-4:
            continue
        # Indices, critical values, and values before and after potential spike
        im1, i0, ip1, ip2, ip3, ip4 = [i-1, i, i+1, i+2, i+3, i+4]
        tm1, t0, tp1, tp2, tp3, tp4 = diff.iloc[[im1, i0, ip1, ip2, ip3, ip4]]
        cm1, c0, cp1, cp2, cp3, cp4 = crit.iloc[[im1, i0, ip1, ip2, ip3, ip4]]
        # print(dates.iloc[[i0,ip1]])
        # print(" sgn(t0)[{}] != sgn(tm1)[{}] \n abs(tm1)[{}] < cm1[{}] \n abs(tp1)[{}] > cp1[{}] \n abs(tp2)[{}] < cp2[{}]".format(np.sign(t0),np.sign(tp1), np.abs(tm1), cm1, np.abs(tp1), cp1, np.abs(tp2), cp2))
        # print()
        # Three-values spike
        if (
            np.sign(t0) != np.sign(tp2) and 
            np.abs(tm1) < 0.5*cm1 and 
            np.abs(tp1) < 0.5*cp1 and 
            np.abs(tp2) < 0.5*cp2 and 
            np.abs(tp3) > cp3 and 
            np.abs(tp4) < 0.5*cp4
        ):
            # print(dates.iloc[i0])
            spikes.iloc[[i0,ip1,ip2]] = True
            # i += 3
            # continue
            
        # Two-values spike
        elif (
            np.sign(t0) != np.sign(tp2) and 
            np.abs(tm1) < 0.5*cm1 and 
            np.abs(tp1) < 0.5*cp1 and 
            np.abs(tp2) > cp2 and 
            np.abs(tp3) < 0.5*cp3
        ):
            # print(dates.iloc[i0])
            spikes.iloc[[i0,ip1]] = True
            # i += 2
            # continue
        
        # One-value spike
        elif( 
            np.sign(t0) != np.sign(tp1) and 
            np.abs(tm1) < 1.0*cm1 and 
            np.abs(tp1) > cp1 and 
            np.abs(tp2) < 1.0*cp2
        ):
            # print(dates.iloc[i0])
            spikes.iloc[i0] = True
            # i += 1
            # continue
        
    return spikes

#---------------------------------------------------------------------------------------------------
def detect_spikes(df, var, iqr_thresh=6, min_datapoints=50):
    """
    Detect  unusual large jumps or ''spikes'' in the time series for `var`:
      1- Find potential unusual large jumps or ''spikes'' by comparing the differences in `var` to each 
         month's critical value (crit = iqr_thresh * IQR)
      2- `potential_spike_check` checks for neccessary conditions for a potential spike to be an actual
          spike
    
    This test is done for ["tas", "tdps", "ps", "slp"]
    Should it be done for more vars?
    
    Input:
    -----
            df [pandas dataframe] : station dataset converted to dataframe through QAQC pipeline
            var [str] : variable to test
            iqr_thresh [int] : critical value (iqr_thresh*IQR) for spike detection (default=6)
            min_datapoints [int] : minimum data points in each month to be valid for testing (default=50)
    Output:
    ------
            df [pandas dataframe] : input df with added columns for spike check
            
            
    NOTES (TODO:)
    iqr_thresh is something can me modified of tweaked down the line (6 is what HadISD uses)
    min_datapoints is the minimum data points in a group for threshold calculation (month/hours between data points)
    HadISD uses 100, this can be modified and twaked in future development
    """
    
    # Make a copy of the original dataframe
    df = df.copy(deep=True)
    
    # Calculate difference in var values
    df[var+'_difference'] = df[var].diff().fillna(0)
    
    # Calculate dates
    df['date'] = df.index.values
    
    # Calculate time difference
    df['time_diff'] = df['date'].diff().fillna(pd.Timedelta(0))
    
    # Calculate time differece in hours
    df['hours_diff'] = df['time_diff']/np.timedelta64(1, 'h')
    df = df[np.logical_and(df['hours_diff'] > 0, df['hours_diff'] <= 12)]
    
    # Group by month to avoid strong seasonal cycle
    # grouped = df.groupby([pd.Grouper(freq='M'), df['hours_diff']]) 
    grouped = df.groupby(pd.Grouper(freq='M')) 

    # Count number of data per month
    counts = grouped[var+'_difference'].transform("count")
    df[var+'_counts'] = counts
    # Keep only months with more than 50 values to be statistically valid
    df = df[df[var+'_counts']>min_datapoints]

    # Define modified IQR 
    # kwargs = {'rng':(20, 80),}
    kwargs = {}
    
    # Calculate iqr
    iqr = grouped[var+'_difference'].transform(stats.iqr, **kwargs)
    df[var+'_iqr'] = iqr

    # Calculate critical value as rounded-up 6 (or defined by argument) times IQR
    df[var+'_critical'] = np.ceil(iqr_thresh*df[var+'_iqr'])
    
    # Find potential spike values where var diff is higher than the critical value
    df[var+'_potential_spikes'] = np.abs(df[var+'_difference'])>df[var+'_critical']

    # Filter real spikes using `potential_spike_check` function
    spikes = potential_spike_check(df[var+'_potential_spikes'], df[var+'_difference'], df[var+'_critical'], df['hours_diff'])
    df[var+'_spikes'] = spikes
    
    return df

#---------------------------------------------------------------------------------------------------
def qaqc_unusual_large_jumps(df, iqr_thresh=6, min_datapoints=50, plot=True, verbose=True):
    """
    Test for unusual large jumps or ''spikes'', given the statistics of the series. Analysis for each individual month in 
    time series to account for seasonal cycles in different regions.
    
    This test is done for ["tas", "tdps", "ps", "psl", "ps_altimeter"]
    Should it be done for more vars?
    
    Input:
    -----
            df [pandas dataframe] : station dataset converted to dataframe through QAQC pipeline
            iqr_thresh [int] : critical value (iqr_thresh*IQR) for spike detection (default=6)
            min_datapoints [int] : minimum data points in each month to be valid for testing (default=50)
            plot [bool] : if True, produces plot and uploads it to AWS
    Output:
    ------
            if qaqc succeded:
                 df [pandas dataframe] : QAQC dataframe with flagged values (see below for flag meaning).
            else if qaqc failed:
                 None
    Flag meaninig:
    -------------
        22,qaqc_unusual_large_jumps,Unusual jump (spike) in variable
    """
    
    try:
        # Save original df multiindex and create station column
        MultiIndex = df.index
        df['station'] = df.index.get_level_values(0)

        # Drop station index
        df = df.droplevel(level="station")

        # Define test variables and check if they are in the dataframe
        check_vars = ["tas", "tdps", "tdps_derived", "ps", "psl", "ps_altimeter"]
        variables = [var for var in check_vars if var in df.columns]

        if verbose:
            print("Running {} on {}".format("qaqc_unusual_large_jumps", variables))

        # Loop through test variables
        for var in variables:
            # Create a copy of the original dataframe and drop NaNs in the testing variable
            new_df = df.copy(deep=True)
            new_df = new_df.dropna(subset=var)#.drop(columns=["lat","lon","elevation"])
            # Use only values that have not been flagged by previous QAQC tests
            valid = np.where(np.isnan(new_df[var+"_eraqc"]))[0]
            new_df = new_df.iloc[valid]
            # Detect spikes
            new_df = detect_spikes(new_df, var=var, iqr_thresh=iqr_thresh, min_datapoints=min_datapoints)
            # Retrieve location of spikes
            ind = new_df.index[np.where(new_df[var+"_spikes"])[0]]
            # Flag _eraqc variable
            df.loc[ind, var+"_eraqc"] = 22

            if plot:
                unusual_jumps_plot(df, var, flagval=22)

        # Reassign original MultiIndex
        df.set_index(MultiIndex, inplace=True)

        return df
    except Exception as e:
        print("qaqc_unusual_large_jumps failed with Exception: {}".format(e))
        return None

# UNUSUAL LARGE JUMPS DEV
#===================================================================================================

#----------------------------------------------------------------------
# To do
# establish false positive rate
# unit tests on each of these functions(?)
