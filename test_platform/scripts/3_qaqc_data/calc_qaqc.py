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
        print('station does not report precipitation - bypassing precip logic nonnegvals check')
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
    
    # pr: Precipitation accumulated since last record
    # pr_5min: Precipitation accumulated in last 5 minutes
    # pr_1h: Precipitation accumulated in last hour
    # pr_24h: Precipitation accumulated from last 24 hours
    # pr_localmid: Precipitation accumulated from local midnight
        
    # rules
    # pr_5min < pr_1h < pr_24h
    # pr_localmid should never exceed pr_24h
    """

    # identify which precipitation vars are reported by a station
    all_pr_vars = [var for var in df.columns if 'pr' in var] # can be variable length depending if there is a raw qc var
    pr_vars = [var for var in all_pr_vars if 'qc' not in var] # remove all qc variables so they do not also run through: raw, eraqc, qaqc_process
    pr_vars = [var for var in pr_vars if 'method' not in var]
    pr_vars = [var for var in pr_vars if 'duration' not in var]

    if not pr_vars: # precipitation variable(s) is not present
        print('station does not report precipitation - bypassing precip logic accum check')
        return None
    
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
    
    # Extract station name
    station = df['station'].unique()[0]
    
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
    df['month'] = pd.to_datetime(df['time']).dt.month # sets month to new variable
    df['year'] = pd.to_datetime(df['time']).dt.year # sets year to new variable
    
    # calculate monthly medians
    df_anom = df.sub(df.resample('M', on='time').transform('median', numeric_only=True))
    df_anom['time'] = df['time'] # add time column back in to do quantiles

    # standardize against calendar-month IQR range
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
#----------------------------------------------------------------------
# To do
# establish false positive rate
# unit tests on each of these functions(?)
