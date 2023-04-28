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

    return dem_elev_short


def qaqc_elev_infill(df):
    """
    Checks if elevation is NA/missing. If missing, fill in elevation from either DEM or station.
    Some stations have all nan elevation values (e.g., NDBC, MARITIME)
    Some stations have single/few but not all nan elevation values (e.g., otherisd, asosawos)
    """

    # create "change" column to check that lat-lon has not changed over time
    df['lat_change'] = df['lat'].shift() != df['lat']
    df['lon_change'] = df['lon'].shift() != df['lon']

    # first check to see if any elev value is missing
    if df['elevation'].isnull().any() == True: 
        print('This station reports a NaN for elevation, infilling from DEM') # testing

        # all elevation values are reported as nan (some ndbc/maritime)
        if df['elevation'].isnull().values.all() == True: 
            try:  # In-fill if value is missing

                # check to see if lat or lon has changed
                if (df['lat_change'][1:].any() == True) or (df['lon_change'][1:].any() == True): # lat or lon has shifted, run through all idxs
                    for i in range(df.shape[0]):
                        dem_elev_value = _grab_dem_elev_m(df['lat'].iloc[i], df['lon'].iloc[i]) # infill dem value
                        df.loc[df.index[i], 'elevation'] = dem_elev_value
                        df.loc[df.index[i], 'elevation_eraqc'] = '3'    ## QC FLAG FOR DEM FILLED VALUE

                else: # lat or lon has not changed, can infill all from iloc[0]
                    dem_elev_value = _grab_dem_elev_m(df['lat'].iloc[0], df['lon'].iloc[0])
                    df['elevation'] = df['elevation'].fillna(dem_elev_value) 
                    df['elevation_qc'] = df["elevation_qc"].fillna('3')   

            except: # elevation cannot be obtained from DEM
                df = pd.DataFrame() # returns empty df to flag that it does not pass


        # some stations have a single/few nan reported (some otherisd/asosawos stations)
        else: 
            print('this station has a missing elevation, but not all -- in progress') # dummy flag message to note which stations this occurs for focused testing

            # if elev changes (is nan) over time, and lat/lon also change (station has moved)
            for i in range(df.shape[0]):
                try:
                    if df['elevation'].iloc[i] != df['elevation'].iloc[0]: # elevation value has changed over time
                        if (df['lat_change'][1:].any() == True) or (df['lon_change'][1:].any() == True): # lat and lon have changed -- pull from dem
                            elev_to_fill = _grab_dem_elev_m(df['lat'].iloc[i], df['lon'].iloc[i])
                            df.loc[df.index[i], 'elevation'] = elev_to_fill
                            df.loc[df.index[i], 'elevation_eraqc'] = '3'
                        
                        else: # elev changes (nan), but lat lon does not -- pull from dem
                            stn_elev = df['elevation'].iloc[0]
                            df.loc[df.index[i], 'elevation'] = stn_elev
                            df.loc[df.index[i], 'elevation_eraqc'] = '4' ## QC FLAG FOR STATION FILLED VALUE
                
                except: # elevation cannot be obtained from DEM, or infilled from station
                    df = pd.DataFrame() # returns empty df to flag that it does not pass

    else:
        df = df

    # drop change columns - not needed further
    df.drop(columns=['lat_change', 'lon_change'], inplace=True)

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


def qaqc_elev_consistency(df):
    """
    Checks if valid elevation value changes drastically from its reported values.
    If drastic change, time record is flagged as suspect.
    """

    # consistency check
    # if elevation value is valid, but significantly different than "normal"
    # e.g. ASOSAWOS_72479723176 elev is 1534m, but has valid records of 0m (flagged in the raw_qc)
    
    # create "change" column to check that elevation has not changed over time
    df['elev_change'] = df['elevation'].shift() != df['elevation']
    base_elev = df['elevation'].iloc[0]

    if (df['elev_change'][1:].any() == True):
        # compare to first value (or value prior)
        for i in range(df.shape[0]):
            if np.abs(df['elevation'][i] - base_elev) > 50: # arbritary threshold of 50 m here --> change to 10/20% different

                # check with dem first
                elev_to_check = _grab_dem_elev_m(df['lat'].iloc[i], df['lon'].iloc[i])

                if np.abs(df['elevation'][i] - float(elev_to_check)) > 5: # arbitrary threshold of 5 m here, but should be close to DEM value
                    continue # elevation does match lat-lon

                else:
                    df.loc[df.index[i], 'elevation_eraqc'] = '1' ## QC FLAG FOR SUSPECT


    ### think about repeating values here
    
    else: # elevation consistent throughout record
        df = df

    # drop change columns - not needed further
    df.drop(columns=['elev_change'], inplace=True)

    return df



                
    ## changes between 0 and 1500m (example), but not nan
        ## check lat lon value from DEM




## Time conversions
## Need function to calculate sub-hourly to hourly -- later on?

#----------------------------------------------------------------------
## Part 2 functions (individual variable/timestamp)


#----------------------------------------------------------------------
# To do
# establish false positive rate
# unit tests on each of these functions(?)
