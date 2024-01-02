"""
This is a script where Stage 3: QA/QC function(s) whole station checks. 
For use within the PIR-19-006 Historical Obsevations Platform.
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
import shapely
import xarray as xr
import matplotlib.pyplot as plt
from io import BytesIO, StringIO
import scipy.stats as stats

wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"

#======================================================================
## Part 1a functions (whole station/network)
## Note: QA/QC functions in part 1a of whole station checks do not proceed through QA/QC if failure occurs

# missing value check: double check that all missing value observations are converted to NA before QA/QC
def qaqc_missing_vals(df, verbose=True):
    '''
    Test for any errant missing values that made it through cleaning and converts missing values to NaNs.
    Searches for missing values in 3_qaqc_data/missing_data_flags.csv. 

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline

    Output:
    -------
        if qaqc success:
            df [pd.DataFrame]: QAQC dataframe with missing values replaced
        if qaqc failure:
            None
    '''

    missing_vals = pd.read_csv('missing_data_flags.csv')

    all_vars = [col for col in df.columns if 'qc' not in col]
    obs_vars = [var for var in all_vars if var not in ['lon','lat','time','elevation','station','anemometer_height_m','thermometer_height_m']]
    
    try:
        for item in obs_vars:
            # pull missing values which are appropriate for the range of real values for each variable 
            missing_codes = missing_vals.loc[missing_vals['variable'].str.contains(item) | missing_vals['variable'].str.contains('all')]

            # values in column that == missing_flag values, replace with NAs
            # note numerical vals converted to strings first to match missing_flag formatting
            df[item] = np.where(df[item].astype(str).isin(missing_codes['missing_flag']), float('NaN'), df[item])

            print('Updating missing values for: {}'.format(item))
    except:
        return None

    return df

#----------------------------------------------------------------------
# missing spatial coords (lat-lon)
def qaqc_missing_latlon(df, verbose=True):
    '''
    Test for missing latitude / longitude values for a station. 
    Checks if latitude and longitude is missing for a station.
    If missing, station is flagged to not proceed through QA/QC.

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline

    Output:
    -------
        if qaqc success:
            df [pd.DataFrame]: QAQC dataframe with missing values replaced
        if qaqc failure:
            None: station does not proceed through QA/QC
    '''

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
    '''
    Test for whether station is within terrestrial & marine WECC boundaries.
    If outside of boundaries, station is flagged to not proceed through QA/QC.

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline

    Output:
    -------
        if qaqc success:
            df [pd.DataFrame]: QAQC dataframe with missing values replaced
        if qaqc failure:
            None: station does not proceed through QA/QC
    '''
    
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
    '''
    If elevation is missing, retrieves elevation value from the USGS Elevation Point Query Service, 
    lat lon must be in decimal degrees (which it is after cleaning). 
    
    Input:
    ------
        lats_to_check [list]: list of latitude values to input for DEM elevation retrieval
        lons_to_check [list]: list of longitude values to input for DEM elevation retrieval

    Output:
    -------
        dem_elev_short [float]: elevation infill value

    References:
    -----------
    Modified from: https://gis.stackexchange.com/questions/338392/getting-elevation-for-multiple-lat-long-coordinates-in-python
    '''
    
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
        dem_elev_short[i] = np.round(dem_elev_long, decimals=2) 

    return dem_elev_short.astype("float")

#----------------------------------------------------------------------
def qaqc_elev_infill(df, verbose=True):
    '''
    Test if elevation is NA/missing. 
    Three in-fill scenarios:
        1 - If all values are nan, in-fill from from DEM
        2 - If some values are missing, in-fill from station for consistency (e.g., OtherISD, ASOSAWOS)
        3 - If all values are nan, and DEM fails to retrieve (over ocean), set to 0 (e.g., NDBC, MARITIME)

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline

    Output:
    -------
        if QAQC success:
            df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        if failure:
            None: station does not proceed through QA/QC

    Flag meaning:
    -------------
        3,qaqc_elev_infill,Elevation infilled from DEM (USGS 3DEP)
        4,qaqc_elev_infill,Elevation infilled from station
        5,qaqc_elev_infill,Elevation manually infilled for buoy out-of-range of DEM

    Note:
    -----
        Elevation qa/qc flags: observations with these flags should continue through QA/QC
    '''
    
    if verbose:
        print('Elevation values pre-infilling: {}'.format(df['elevation'].unique()))
        print('Elevation eraqc values pre-infilling: {}'.format(df['elevation_eraqc'].unique())) # testing

    # first check to see if any elev value is missing
    if df['elevation'].isnull().any() == True: 
        # all elevation values are reported as nan (some ndbc/maritime)

        if df['elevation'].isnull().values.all() == True: 
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
                        df.loc[df['elevation'].isnull(), 'elevation'] = float(0.00) # manual infilling for buoys

                else: # multiple pairs of lat-lon for missing elevs
                    for ilat in nan_lats:
                        for ilon in nan_lons:
                            dem_elev_value = _grab_dem_elev_m(ilat, ilon)
                            df.loc[(df['lat'] == ilat) & (df['lon'] == ilon), 'elevation_eraqc'] = 3 # see era_qaqc_flag_meanings.csv
                            df.loc[(df['lat'] == ilat) & (df['lon'] == ilon), 'elevation'] = float(dem_elev_value)
                        
            except: # elevation cannot be obtained from DEM
                if verbose:
                    print("Elevation cannot be in-filled")
                return None

        # some stations have a single/few nan reported, in-fill from station for consistency
        else:   # multiple values for elevation, infill each instance if missing/incorrectly coded (e.g., zeros when shouldnt be)
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
                                                            
            # elevation cannot be in-filled
            except:
                if verbose:
                    print("Elevation cannot be in-filled")
                return None
    else:
        return df

#----------------------------------------------------------------------
def qaqc_elev_range(df, verbose=True):
    '''
    Checks if valid elevation value is outside of range of reasonable values for WECC region.
    If outside range, station is flagged to not proceed through QA/QC.

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline

    Output:
    -------
        if QAQC success:
            df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        if failure:
            None
    '''

    # death valley is 282 feet (85.9 m) below sea level
    # denali is ~6190 m
    # add a 10 m buffer

    # If value is present but outside of reasonable value range
    if (df['elevation'].values.any() < -95.0) or (df['elevation'].values.any() > 6210.0):
        if verbose:
            print("Station out of range for WECC -- station does not proceed through QAQC")
            return None

    # Elevation value is present and within reasonable value range
    else:
        df = df
        if verbose:
            print('Elevation values post-infilling/correcting: {}'.format(df['elevation'].unique())) # testing
            print('Elevation qaqc values post-infilling/correcting: {}'.format(df['elevation_eraqc'].unique())) # testing
        
    return df

#----------------------------------------------------------------------
## sensor height - air temperature
def qaqc_sensor_height_t(df, verbose=True):
    '''
    Checks if temperature sensor height is within 2 meters above surface +/- 1/3 meter tolerance.
    If missing or outside range, temperature value for station is flagged to not proceed through QA/QC.

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline

    Output:
    -------
        if QAQC success:
            df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        if failure:
            None

    Flag meaning:
    -------------
        6,qaqc_sensor_height_t,Thermometer height missing
        7,qaqc_sensor_height_t,Thermometer height not 2 meters
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
                # df.loc[:, 'tas_eraqc'] = 7 # see era_qaqc_flag_meanings.csv
                df['tas_eraqc'] = 7 # see era_qaqc_flag_meanings.csv

        return df
    except Exception as e:
        if verbose:
            print("qaqc_sensor_height_w failed with Exception: {}".format(e))
        return None

#----------------------------------------------------------------------
## sensor height - wind
def qaqc_sensor_height_w(df, verbose=True):
    '''
    Checks if wind sensor height is within 10 meters above surface +/- 1/3 meter tolerance.
    If missing or outside range, wind speed and direction values for station are flagged to not proceed through QA/QC.

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline

    Output:
    -------
        if QAQC success:
            df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        if failure:
            None

    Flag meaning:
    -------------
        8,qaqc_sensor_height_w,Anemometer height missing
        9,qaqc_sensor_height_w,Anemometer height not 10 meters
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

#======================================================================
## Part 1b functions (whole station/network)
## Note: QA/QC functions in part 1b of whole station checks proceed through QA/QC if failure occurs

## flag values outside world records for North America
def qaqc_world_record(df, verbose=True):
    '''
    Checks if temperature, dewpoint, windspeed, or sea level pressure are outside North American world records
    If outside minimum or maximum records, flags values. 

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline

    Output:
    -------
        if QAQC success:
            df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        if failure:
            None

    Flag meaning:
    -------------
        11,qaqc_world_record,Value outside of world record range
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
                    df.loc[isOffRecord, var + '_eraqc'] = 11 # see era_qaqc_flag_meanings.csv
        return df
    except Exception as e:
        if verbose:
            print("qaqc_world_record failed with Exception: {}".format(e))
        return None
