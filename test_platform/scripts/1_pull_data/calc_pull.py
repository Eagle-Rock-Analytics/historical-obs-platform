"""
This is a script where Stage 1: Pull related common functions, conversions, and operations is stored for ease of use
for the Historical Observations Platform.
"""

## Process to call this script: from SCRIPT import FUNCTION
## Example: from calc_pull import get_wecc_poly

## Import Libraries
import geopandas as gp


## Useful functions
def get_wecc_poly(terrpath, marpath):
    """
    Identifies a bbox of WECC area to filter stations against
    Input vars: shapefiles for maritime and terrestrial WECC boundaries
    Returns: spatial objects for each shapefile, and bounding box for their union.
    """
    t = gp.read_file(terrpath)  ## Read in terrestrial WECC shapefile.
    m = gp.read_file(marpath)  ## Read in marine WECC shapefile.
    bbox = t.union(m).bounds  ## Combine polygons and get bounding box of union.
    return t, m, bbox


def _lat_dms_to_dd(data):
    """
    Converts latitude from decimal-minutes-seconds to decimal degrees
    Input: latitude (DMS) example: 39.02.33
    Returns: latitude (dd) example: 39.16
    """
    data = float(data[:2]) + float(data[3:5]) / 60 + float(data[6:]) / 3600
    return data


def _lon_dms_to_dd(data):
    """
    Converts longitude from decimal-minutes-seconds to decimal degrees
    and ensures that western hemisphere lons are negative by convention
    Input: longitude(DMS) example: 122.01.38
    Returns: longitude (dd) example: -122.02
    """
    # need to check if -180 to 180, or 0 to 360
    if data[0] != "-":
        _deg = float(data[:3])
        _min = float(data[4:6])
        _sec = float(data[7:])
        data = -1 * (_deg + _min / 60 + _sec / 3600)
    else:
        data = data.strip("-")
        _deg = float(data[:3])
        _min = float(data[4:6])
        _sec = float(data[7:])
        data = -1 * (_deg + _min / 60 + _sec / 3600)
    return data
