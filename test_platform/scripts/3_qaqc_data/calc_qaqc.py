"""
This is a script where Stage 3: QA/QC related common functions, conversions, and operations is stored for ease of use
for the Historical Observations Platform.
"""

## Import Libraries
import boto3
import geopandas as gp
import numpy as np
import pandas as pd


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
def qaqc_elev_demfill(file_to_qaqc):
    """
    Checks if elevation is NA/missing. If missing, fill in elevation from DEM.
    """

    # elevation is nan
    if file_to_qaqc['elevation'].isnull().any() == True: 
        if file_to_qaqc['elevation'].isnull().values.all() == True: # all elevation values are reported as nan (some ndbc/maritime)
            file_to_qaqc['elevation_qc'] = file_to_qaqc["elevation_qc"].fillna("E")   ## FLAG FOR DEM FILLED VALUE
            # try:  # In-fill if value is missing
                # dem = rio.open(dem)
                # dem_array = dem.read(1).astype('float64')
                # make sure to round off lat-lon values so they are not improbably precise for our needs

            # except:
                # file_to_qaqc = pd.DataFrame() # returns empty df to flag that it does not pass
            print('This station reports a NaN for elevation, infilling from DEM -- in progress')
            file_to_qaqc = file_to_qaqc
            # file_to_qaqc = pd.DataFrame() # returns empty df to flag that it does not pass

        # else: 
        #     print('test') # some stations have a single nan reported (some otherisd)
        #     identify which rows/timestamps have a missing elevation value
        #     infill from dem

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
