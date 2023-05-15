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
                        df.loc[df['elevation'].isnull(), 'elevation_eraqc'] = 3 # see qaqc_flag_meanings.csv
                        df.loc[df['elevation'].isnull(), 'elevation'] = float(dem_elev_value)
                    except: # some buoys out of range of dem (past coastal range) report nan elevation, manually set to 0.00m and flag
                        df.loc[df['elevation'].isnull(), 'elevation_eraqc'] = 5 # see qaqc_flag_meanings.csv
                        df.loc[df['elevation'].isnull(), 'elevation'] = float(0.00) # manual infilling

                else: # multiple pairs of lat-lon for missing elevs
                    for ilat in nan_lats:
                        for ilon in nan_lons:
                            dem_elev_value = _grab_dem_elev_m(ilat, ilon)
                            df.loc[(df['lat'] == ilat) & (df['lon'] == ilon), 'elevation_eraqc'] = 3 # see qaqc_flag_meanings.csv
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
                        df.loc[df['elevation'].isnull(), 'elevation_eraqc'] = 4 # see qaqc_flag_meanings.csv
                        df.loc[df['elevation'].isnull(), 'elevation'] = df['elevation'].iloc[0]

                    else: # lat-lon of missing elev does not match station lat-lon (has shifted), infill from dem
                        dem_elev_value = _grab_dem_elev_m(nan_lats[0], nan_lons[0])
                        df.loc[df['elevation'].isnull(), 'elevation_eraqc'] = 3 # see qaqc_flag_meanings.csv
                        df.loc[df['elevation'].isnull(), 'elevation'] = float(dem_elev_value)

                else: # multiple pairs of lat-lon for missing elevs
                    for ilat in nan_lats:
                        for ilon in nan_lons:
                            dem_elev_value = _grab_dem_elev_m(ilat, ilon)
                            df.loc[(df['lat'] == ilat) & (df['lon'] == ilon), 'elevation_eraqc'] = 3 # see qaqc_flag_meanings.csv
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
                    df.loc[df['elevation'] == elev_value, 'elevation_eraqc'] = 4 # see qaqc_flag_meanings.csv
                    df.loc[df['elevation'] == elev_value, 'elevation'] = df['elevation'].iloc[0]

                else: # lat-lon of incorrectly coded elevation does not match station lat-lon (has shifted), infill from dem
                    dem_elev_value = _grab_dem_elev_m(off_lats[0], off_lons[0])
                    df.loc[df['elevation'] == elev_value, 'elevation_eraqc'] = 3 # see qaqc_flag_meanings.csv
                    df.loc[df['elevation'] == elev_value, 'elevation'] = float(dem_elev_value)
            
            else: # multple pairs of lat-lon for incorrectly zero coded elevs 
                for ilat in off_lats:
                    for ilon in off_lons:
                        dem_elev_value = _grab_dem_elev_m(ilat, ilon)
                        df.loc[(df['lat'] == ilat) & (df['lon'] == ilon), 'elevation_eraqc'] = 3 # see qaqc_flag_meanings.csv
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
    pr_vars = [var for var in all_pr_vars if "_qc" not in var] # remove all qc variables so they do not also run through

    if len(pr_vars) != 0: # precipitation variable(s) is present
        for item in pr_vars:
            if (df[item] < 0).any():
                df.loc[df[item] < 0, 'pr_eraqc'] = 10 # see qaqc_flag_meanings.csv

    else: # station does not report precipitation
        print('station does not report precipitation - bypassing precip logic check') # testing
        df = df

    print('Precipitation neg-values logic check testing: {}'.format(df['pr_eraqc'].unique())) # testing

    return df




#----------------------------------------------------------------------
# To do
# establish false positive rate
# unit tests on each of these functions(?)
