"""
This is a script where Stage 1: Pull related common functions, conversions, and operations is stored for ease of use
for the Historical Observations Platform.
"""

## Process to call this script: from SCRIPT import FUNCTION
## Example: from calc_pull import get_wecc_poly

## Import Libraries
import numpy as np
import geopandas as gp
from math import exp, log


## Useful functions
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
