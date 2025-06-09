"""
calc_pull.py

This is a script where Stage 1: Pull related common functions, conversions, and operations is stored for ease of use
for the Historical Observations Platform.

Functions
---------
- get_wecc_poly: Identifies a bbox of WECC area to filter stations against
- _lat_dms_to_dd: Converts latitude from decimal-minutes-seconds to decimal degrees
- _lon_dms_to_dd: Converts longitude from decimal-minutes-seconds to decimal degrees
"""

import geopandas as gpd
from shapely.geometry import box

def get_wecc_poly(terrpath: str, marpath: str) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoSeries]:
    """
    Identifies a bbox of WECC area to filter stations against

    Parameters
    ----------
    terrpath : str
        shapefiles for maritime and terrestrial WECC boundaries
    marpath : str
        shapefiles for maritime and terrestrial WECC boundaries

    Returns
    -------
    t : gpd.GeoDataFrame
        spatial object of terrestrial WECC
    m : gpd.GeoDataFrame
        spatial object of maritime WECC
    bbox : gpd.GeoSeries
        spatial object bounding box
    """
    t = gpd.read_file(terrpath) 
    m = gpd.read_file(marpath) 

    # Combine polygons and get bounding box of union
    combined = t.geometry.unary_union.union(m.geometry.unary_union)
    bbox = gpd.GeoSeries([box(*combined.bounds)], crs=t.crs)  

    return t, m, bbox


def _lat_dms_to_dd(data: str) -> float:
    """
    Converts latitude from decimal-minutes-seconds to decimal degrees

    Parameters
    ----------
    data : str
        input latitude

    Returns
    -------
    data : float
        converted latitude

    Example
    -------
    latitude (DMS) input: 39.02.33
    latitude (dd) output: 39.16
    """
    data = float(data[:2]) + float(data[3:5]) / 60 + float(data[6:]) / 3600
    return data


def _lon_dms_to_dd(data: str) -> float:
    """
    Converts longitude from decimal-minutes-seconds to decimal degrees
    and ensures that western hemisphere lons are negative by convention

    Parameters
    ----------
    data : str
        input longitude

    Returns
    -------
    data : float
        converted longitude

    Example
    -------
    longitude(DMS) output: 122.01.38
    longitude (dd) output: -122.02
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
