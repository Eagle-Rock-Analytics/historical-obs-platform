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
def _grab_dem_elev_m(file_to_qaqc):
    """
    Modified from: https://gis.stackexchange.com/questions/338392/getting-elevation-for-multiple-lat-long-coordinates-in-python
    """
    # for USGS EPQS, lat lon must be in decimal degrees (which it is after cleaning)
    lat_to_check = file_to_qaqc['lat'].iloc[0]
    lon_to_check = file_to_qaqc['lon'].iloc[0]

    # USGS Elevation Point Query Service
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


def qaqc_elev_demfill(file_to_qaqc):
    """
    Checks if elevation is NA/missing. If missing, fill in elevation from DEM.
    Some stations have all nan elevation values (e.g., NDBC, MARITIME)
    Some stations have single/few but not all nan elevation values (e.g., otherisd)
    """

    # first check to see if any elev value is missing
    if file_to_qaqc['elevation'].isnull().any() == True: 
        print('This station reports a NaN for elevation, infilling from DEM')

        if file_to_qaqc['elevation'].isnull().values.all() == True: # all elevation values are reported as nan (some ndbc/maritime)
            file_to_qaqc['elevation_qc'] = file_to_qaqc["elevation_qc"].fillna("2")   ## QC FLAG FOR DEM FILLED VALUE
            try:  # In-fill if value is missing
               dem_elev_value = _grab_dem_elev_m(file_to_qaqc)
               file_to_qaqc['elevation'] = file_to_qaqc['elevation'].fillna(dem_elev_value) # infill dem value

            except: # elevation cannot be obtained from DEM
                file_to_qaqc = pd.DataFrame() # returns empty df to flag that it does not pass

        else: # some stations have a single/few nan reported (some otherisd stations) -- come back to in other network testing (use clean_master_station_list!)
            print('this station has a missing elevation, but not all -- in progress') # dummy flag message to note which stations this occurs for focused testing
            dem_elev_value = _grab_dem_elev_m(file_to_qaqc)
            print(dem_elev_value) # testing
            try:
                ## THIS NEEDS TO BE TESTED AND CONFIRMED
                for i in range(file_to_qaqc.shape[0]):       
                    if file_to_qaqc['elevation'] == np.nan: #  identify which rows/timestamps have a missing elevation value
                        file_to_qaqc.loc[file_to_qaqc.index[i], 'elevation'] = dem_elev_value

                        ## notes:
                        ## need to assess what the values before and after the infilled value are and confirm that they match + sig figs
                        ## can alternatively, grab those values instead for consistency, rather than infilling
            
            except:
                file_to_qaqc = pd.DataFrame() # returns empty df to flag that it does not pass

    else:
        file_to_qaqc = file_to_qaqc

    return file_to_qaqc


def qaqc_elev_range(file_to_qaqc):
    """
    Checks if valid elevation value is outside of range of reasonable values for WECC region.
    If outside range, station is flagged to not proceed through QA/QC.
    """

    # death valley is 282 feet (85.9 m) below sea level
    # denali is ~6190 m

    # If value is present but outside of reasonable value range
    if (file_to_qaqc['elevation'].values.any() < -86.0) or (file_to_qaqc['elevation'].values.any() > 6200.0):
        file_to_qaqc = pd.DataFrame() # returns empty df to flag that it does not pass

    # Elevation value is present and within reasonable value range
    else:
        file_to_qaqc = file_to_qaqc

    return file_to_qaqc




## Time conversions
## Need function to calculate sub-hourly to hourly -- later on?

#----------------------------------------------------------------------
## Part 2 functions (individual variable/timestamp)


#----------------------------------------------------------------------
# To do
# establish false positive rate
# unit tests on each of these functions(?)
