"""
This is a script where Stage 2: Clean related common functions, conversions, and operations is stored for ease of use
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


## Unit conversions, if required
## Temperature conversions: Desired working unit should be K
def _unit_degC_to_K(data):
    """
    Converts temperature from degC to K
    Inputs: temperature (degC)
    Returns: temperature (K)
    """
    data = data + 273.15
    return data

def _unit_degF_to_K(data):
    """
    Converts temperature from degF to K
    Inputs: temperature (degF)
    Returns: temperature (K)
    """
    data = (5/9) * (data - 32) + 273.15
    return data

## Air pressure conversions: Desired working unit should be mb or hPa
def _unit_pres_pa_to_hpa(data):
    """
    Converts air pressure from pascals to hectopascals
    Inputs: air pressure (Pa)
    Returns: air pressure (hPa)
    Note: this also works for the conversion to mb
    """
    data = data / 100.
    return data

def _unit_pres_inHg_to_hpa(data):
    """
    Converts air pressure from inHg to hectopascals
    Inputs: air pressure (inHg)
    Returns: air pressure (hPa)
    """
    data = data * 33.8639
    return data

## Wind speed conversions: Desired working unit should be m/s
def _unit_windspd_kts_to_ms(data):
    """
    Converts windspeed from knots to m/s
    Inputs: windspeed (knots)
    Returns: windspeed (m/s)
    """
    data = data / 1.94
    return data

def _unit_windspd_mph_to_ms(data):
    """
    Converts windspeed from miles-per-hour to m/s
    Inputs: windspeed (mph)
    Returns: windspeed (m/s)
    """
    data = data / 2.237
    return data

## Moisture conversions: Desired working unit should be  kg/kg
def _unit_moisture_gkg_to_kgkg(data):
    """
    Converts moisture ratios from g/kg to kg/kg
    Inputs: moisture ratio (g/kg)
    Returns: moisture ratio (kg/kg)
    """
    data = data / 1000
    return data

## Precipitation conversions: Desired working unit should be mm/TIME
## Note on desired unit - the conversion to hourly as a rate will occur in calc_qaqc
def _unit_precip_in_to_mm(data):
    """
    Converts precipitation from inches to mm.
    Inputs: precipitation (inches)
    Returns: precipitation (mm)
    """
    data = data * 25.4
    return data

## Elevation conversions: Desired working unit should be meters
def _unit_elev_ft_to_m(data):
    """
    Converts elevation from feet to meters
    Input: elevation (feet)
    Returns: elevation (m)
    """
    data = data / 0.3048
    return data

##---------------------------------------------------------------------------------------------
## Derived variable calculations
def _calc_dewpointtemp_opt1(tas, hurs):
    """
    Calculates dew point temperature, method 1
    Inputs: air temperature (K), relative humidity (0-100 scale)
    Returns: dew point temperature (K)
    """
    es = 0.611 * exp(5423 * ((1/273) - (1/tas)))   # calculates saturation vapor pressure
    e = (es * hurs)/100                        # calculates vapor pressure, IF NOT ALREADY OBSERVED -- will need ifelse statement
    tdps = ((1/273) - 0.0001844 * log(e/0.611))^-1   # calculates dew point temperature, units = K
    return tdps

def _calc_dewpointtemp_opt2(e):
    """
    Calculates dew point temperature, method 2
    Inputs: vapor pressure (Pa)
    Returns: dew point temperature (K)
    """
    tdps = ((1/273) - 0.0001844 * log(e/0.611))^-1   # calculates dew point temperature, units = K

def _calc_relhumid(tas, tdps):
    """
    Calculate relative humidity
    Inputs: air temperature, dewpoint temperature (both K)
    Returns: realtive humidity (%, or 0-100 scale)
    """
    es = 0.611 * exp(5423 * ((1/273) - (1/tas)))   # calculates saturation vapor pressure using air temp
    e = 0.611 * exp(5423 * ((1/273) - (1/tdps)))   # calculates vapor pressure using dew point temp
    hurs = 100 * (e/es)
    return hurs

def _calc_windmag(u10, v10):
    """
    Calculates wind speed
    Inputs: u and v wind components (both m/s)
    Returns: wind speed/magnitude (m/s)
    """
    sfcWind = np.sqrt((u10)^2  + (v10)^2)   # calculates wind magnitude, units = ms-1
    return sfcWind

def _calc_winddir(u10, v10):
    """
    Calculates wind direction
    Inputs: u and v wind components (both m/s)
    Returns: wind direction (degrees)
    NOTE: this is not functional at the moment due to performance issues
    """
    pass        # this is a complicated calculation -- looking for options
    #return sfcWind_dir

## CHECK THIS CALCULATION
def _calc_ps(psl, elev, temp):
    """
    Calculates station air pressure from sea level air pressure, if station pressure is not available
    Inputs: sea level pressure (mb/hPa), elevation (m), and air temperature (K)
    Returns: air pressure (mb/hPa)
    Note: this calculation checks out using 2 different formula, with differences at the second decimal place:
    https://keisan.casio.com/exec/system/1224575267
    https://www.mide.com/air-pressure-at-altitude-calculator
    """
    ps = psl * math.e**(-elev/(temp*29.263))
    return ps
