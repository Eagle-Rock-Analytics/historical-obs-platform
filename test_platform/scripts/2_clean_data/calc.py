'''
This file contains the calculation functions used across all four stages
of the cleaning process for the historical data platform. It should be imported
by other scripts.
'''

# Libraries
import numpy as np
import geopandas as gp
from numpy import sqrt, exp, log

# Spatial operations

## Input vars: shapefiles for maritime and terrestrial WECC boundaries
## Outputs: spatial objects for each shapefile, and bounding box for their union.
def get_wecc_poly(terrpath, marpath):
    ## get bbox of WECC to use to filter stations against
    ## Read in terrestrial WECC shapefile.
    t = gp.read_file(terrpath)
    ## Read in marine WECC shapefile.
    m = gp.read_file(marpath)
    ## Combine polygons and get bounding box of union.
    bbox = t.union(m).bounds
    return t,m, bbox


# Derived variable calculations

## dew point temperature calculation 
## (necessary input vars: requires at least 2 of three - air temp + relative humidity + vapor pressure)
def _calc_dewpointtemp(tas, hurs, e):
    es = 0.611 * exp(5423 * ((1/273) - (1/tas)))   # calculates saturation vapor pressure
    e = (es * hurs)/100                        # calculates vapor pressure, IF NOT ALREADY OBSERVED -- will need ifelse statement
    tdps = ((1/273) - 0.0001844 * log(e/0.611))^-1   # calculates dew point temperature, units = K
    return tdps

# relative humidity calculation (necessary input vars: air temp + dew point**, air temp + vapor pressure, air pressure + vapor pressure)
def _calc_relhumid(tas, tdps):
    es = 0.611 * exp(5423 * ((1/273) - (1/tas)))   # calculates saturation vapor pressure using air temp
    e = 0.611 * exp(5423 * ((1/273) - (1/tdps)))   # calculates vapor pressure using dew point temp
    hurs = 100 * (e/es)
    return hurs

# wind speed (necessary input vars: u and v components)
def _calc_windmag(u10, v10):
    sfcWind = sqrt((u10)^2  + (v10)^2)   # calculates wind magnitude, units = ms-1
    return sfcWind

# wind direction (necessary input vars: u and v components)
def _calc_winddir(u10, v10):
    pass        # this is a complicated calculation -- looking for options
    #return sfcWind_dir