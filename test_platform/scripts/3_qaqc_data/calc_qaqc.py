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
def qaqc_missing_latlon(file_to_qaqc):
    """
    Checks if latitude and longitude is missing for a station.
    If missing, station is flagged to not proceed through QA/QC.
    """

    # latitude
    if file_to_qaqc['lat'].isnull().values.any() == True:
        file_to_qaqc = None # returns empty df to flag that it does not pass
    else:
        file_to_qaqc = file_to_qaqc

    # longitude
    if file_to_qaqc['lon'].isnull().values.any() == True:
        file_to_qaqc = None # returns empty df to flag that it does not pass
    else:
        file_to_qaqc = file_to_qaqc

    return file_to_qaqc

## in bounds of WECC
def qaqc_within_wecc(file_to_qaqc):
    """
    Checks if station is within terrestrial & marine WECC boundaries.
    If outside of boundaries, station is flagged to not proceed through QA/QC.
    """

    t, m, bbox = get_wecc_poly(wecc_terr, wecc_mar) # Call get_wecc_poly
    lat_to_check = file_to_qaqc['lat'].iloc[0]
    lon_to_check = file_to_qaqc['lon'].iloc[0]

    # latitude
    if (lat_to_check < bbox.miny.values) or (lat_to_check > bbox.maxy.values):
        file_to_qaqc = pd.DataFrame() # returns empty df to flag that it does not pass
    else:
        file_to_qaqc = file_to_qaqc

    # longitude
    if (lon_to_check > bbox.maxx.values) or (lon_to_check < bbox.minx.values):
        file_to_qaqc = pd.DataFrame() # returns empty df to flag that it does not pass
    else:
        file_to_qaqc = file_to_qaqc

    return file_to_qaqc

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

    return float(dem_elev_short)


def qaqc_elev_infill(df):
    """
    Checks if elevation is NA/missing. If missing, fill in elevation from either DEM or station.
    Some stations have all nan elevation values (e.g., NDBC, MARITIME)
    Some stations have single/few but not all nan elevation values (e.g., otherisd, asosawos)
    """

    # first check to see if any elev value is missing
    if df['elevation'].isnull().any() == True: 
        print('This station reports a NaN for elevation, infilling from DEM') # testing

        # all elevation values are reported as nan (some ndbc/maritime)
        if df['elevation'].isnull().values.all() == True: 
            try:  # in-fill if value is missing
                # locate all instances of nan values as elevation codes
                nan_lats = df['lat'].unique()
                nan_lons = df['lon'].unique()

                # check if lat-lon has changed over time
                if (len(nan_lats) == 1) and (len(nan_lons) == 1): # single lat-lon pair for missing elevs
                    dem_elev_value = _grab_dem_elev_m(df['lat'].iloc[0], df['lon'].iloc[0])
                    df.loc[df['elevation'].isnull(), 'elevation_eraqc'] = '3'   # QC FLAG FOR DEM FILLED VALUE
                    df.loc[df['elevation'].isnull(), 'elevation'] = dem_elev_value

                else: # multiple pairs of lat-lon for missing elevs
                    for ilat in nan_lats:
                        for ilon in nan_lons:
                            dem_elev_value = _grab_dem_elev_m(ilat, ilon)
                            df.loc[(df['lat'] == ilat) & (df['lon'] == ilon), 'elevation_eraqc'] = '3'  # QC FLAG FOR DEM FILLED VALUE
                            df.loc[(df['lat'] == ilat) & (df['lon'] == ilon), 'elevation'] = dem_elev_value
                
                df = df
                        
            except: # elevation cannot be obtained from DEM
                df = pd.DataFrame() # returns empty df to flag that it does not pass


        # some stations have a single/few nan reported (some otherisd/asosawos stations)
        else:   # multiple values for elevation, infill each instance if missing/incorrectly coded (zeros when shouldnt be)
            print('This station has a missing elevation, but not all -- in progress') # dummy flag message to note which stations this occurs for focused testing

            elev_vals = df['elevation'].unique() # identify how many different elevation "values" are present

            try:
                # locate all instances of nan values as elevation codes
                nan_coded = df[df['elevation'].isnull()]
                nan_lats = nan_coded['lat'].unique()
                nan_lons = nan_coded['lon'].unique()

                if (len(nan_lats) == 1) and (len(nan_lons) == 1): # single lat-lon pair for missing elevs
                    if (nan_lats[0] == df['lat'].iloc[0]) and (nan_lons[0] == df['lon'].iloc[0]): # single set of lat-lons matches station, can infill for consistency
                        df.loc[df['elevation'].isnull(), 'elevation_eraqc'] = '4' # QC FLAG FOR STATION FILLED VALUE
                        df.loc[df['elevation'].isnull(), 'elevation'] = df['elevation'].iloc[0]

                    else:
                        dem_elev_value = _grab_dem_elev_m(nan_lats[0], nan_lons[0])
                        df.loc[df['elevation'].isnull(), 'elevation_eraqc'] = '3'   # QC FLAG FOR DEM FILLED VALUE
                        df.loc[df['elevation'].isnull(), 'elevation'] = dem_elev_value

                else: # multiple pairs of lat-lon for missing elevs
                    for ilat in nan_lats:
                        for ilon in nan_lons:
                            dem_elev_value = _grab_dem_elev_m(ilat, ilon)
                            df.loc[(df['lat'] == ilat) & (df['lon'] == ilon), 'elevation_eraqc'] = '3'  # QC FLAG FOR DEM FILLED VALUE
                            df.loc[(df['lat'] == ilat) & (df['lon'] == ilon), 'elevation'] = dem_elev_value
                                
                # incorrectly coded as zero elevation when station elevation is not zero
                if (df['elevation'].iloc[0] != 0) and (0 in elev_vals):
                            
                    # locate all instances of incorrect elevation 0 codes
                    zero_coded = df[df['elevation']==0]
                    zero_lats = zero_coded['lat'].unique()
                    zero_lons = zero_coded['lon'].unique()
                    
                    if (len(zero_lats) == 1) and (len(zero_lons) == 1): # single lat-lon pair for bad coded elevs
                        if (zero_lats[0] == df['lat'].iloc[0]) and (zero_lons[0] == df['lon'].iloc[0]): # single set of lat-lons matches station, can infill for consistency
                            df.loc[df['elevation'] == 0, 'elevation_eraqc'] = '4' # QC FLAG FOR STATION FILLED VALUE
                            df.loc[df['elevation'] == 0, 'elevation'] = df['elevation'].iloc[0]
                        else:
                            dem_elev_value = _grab_dem_elev_m(zero_lats[0], zero_lons[0])
                            df.loc[df['elevation'] == 0, 'elevation_eraqc'] = '3' # QC FLAG FOR DEM FILLED VALUE
                            df.loc[df['elevation'] == 0, 'elevation'] = dem_elev_value # is this a mix of types (float vs string)??
                        
                    else: # multple pairs of lat-lon for incorrectly zero coded elevs
                        for ilat in zero_lats:
                            for ilon in zero_lons:
                                dem_elev_value = _grab_dem_elev_m(ilat, ilon)
                                df.loc[(df['lat'] == ilat) & (df['lon'] == ilon), 'elevation_eraqc'] = '3'
                                df.loc[(df['lat'] == ilat) & (df['lon'] == ilon), 'elevation'] = dem_elev_value


            except: # elevation cannot be obtained from DEM
                df = pd.DataFrame() # returns empty df to flag that it does not pass

    else:
        df = df

    return df


def qaqc_elev_range(df):
    """
    Checks if valid elevation value is outside of range of reasonable values for WECC region.
    If outside range, station is flagged to not proceed through QA/QC.
    """

    # death valley is 282 feet (85.9 m) below sea level
    # denali is ~6190 m

    # If value is present but outside of reasonable value range
    if (df['elevation'].values.any() < -86.0) or (df['elevation'].values.any() > 6200.0):
        df = pd.DataFrame() # returns empty df to flag that it does not pass

    # Elevation value is present and within reasonable value range
    else:
        df = df

    return df




## Time conversions
## Need function to calculate sub-hourly to hourly -- later on?

#----------------------------------------------------------------------
## Part 2 functions (individual variable/timestamp)


#----------------------------------------------------------------------
# To do
# establish false positive rate
# unit tests on each of these functions(?)
