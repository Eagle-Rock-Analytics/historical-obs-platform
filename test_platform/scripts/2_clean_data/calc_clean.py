"""
This is a script where Stage 2: Clean related common functions, conversions, and operations is stored for ease of use
for the Historical Observations Platform.
"""

## Process to call this script: from SCRIPT import FUNCTION
## Example: from calc_pull import get_wecc_poly

## Import Libraries
import geopandas as gp
import numpy as np

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

## Air pressure conversions: Desired working unit should be Pa
def _unit_pres_hpa_to_pa(data):
    """
    Converts air pressure from hectopascals to pascals
    Inputs: air pressure (hPa)
    Returns: air pressure (Pa)
    Note: this also works for the conversion from mb
    """
    data = data * 100.
    return data

def _unit_pres_inHg_to_pa(data):
    """
    Converts air pressure from inHg to hectopascals
    Inputs: air pressure (inHg)
    Returns: air pressure (Pa)
    """
    data = data * 3386.39
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
    data = data / 1000.
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
    data = data * 0.3048
    return data

## Latitude/Longitude conversions: Desired working unit should be decimal degrees N/W
## Need to also accomodate inputs such as: 41.56.54
def _lat_dms_to_dd(data):
    """
    Converts latitude from decimal-minutes-seconds to decimal degrees
    Input: latitude (DMS) example: 39.02.33
    Returns: latitude (dd) example: 39.16
    """
    data = float(data[:2]) + float(data[3:5])/60 + float(data[6:])/3600
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
        data = -1 * (_deg + _min/60 + _sec/3600)
    else:
        data = data.strip('-')
        _deg = float(data[:3])
        _min = float(data[4:6])
        _sec = float(data[7:])
        data = -1 * (_deg + _min/60 + _sec/3600)
    return data

def _lon_DMm_to_Dd(data):
    """
    This is specific to CWOP longitude data converting from LORAN (DM.m) coordinates to decimal-degrees (D.d) for the WESTERN HEMISPHERE.
    Input: longitude (DDDMM.mm) example: 12234.72
    Returns: longitude (D.d) example: -122.578
    """
    _min = float(data[:3])
    _sec = float(data[3:])
    data = -1 * (_min + _sec/60)
    return data

def _lat_DMm_to_Dd(data):
    """
    This is specific to CWOP latitude data converting from LORAN (DM.m) coordinates to decimal-degrees (D.d).
    Input: latitude (DDMM.mm) example: 4413.95
    Returns: latitude (D.d) example: 44.2325
    """
    _deg = float(data[:2])
    _mm = float(data[2:])
    data = (_deg + _mm/60)
    return data

##---------------------------------------------------------------------------------------------
## Derived variable calculations
def _calc_dewpointtemp_opt1(tas, hurs):
    """
    Calculates dew point temperature, method 1
    Inputs: air temperature (K), relative humidity (0-100 scale)
    Returns: dew point temperature (K)
    """
    es = 0.611 * np.exp(5423 * ((1/273) - (1/tas)))   # calculates saturation vapor pressure
    e_vap = (es * hurs)/100.                        # calculates vapor pressure, IF NOT ALREADY OBSERVED -- will need ifelse statement
    tdps = ((1/273) - 0.0001844 * np.log(e_vap/0.611))**-1   # calculates dew point temperature, units = K
    return tdps

def _calc_dewpointtemp_opt2(e_vap):
    """
    Calculates dew point temperature, method 2
    Inputs: vapor pressure (Pa)
    Returns: dew point temperature (K)
    """
    tdps = ((1/273) - 0.0001844 * np.log(e_vap/0.611))**-1   # calculates dew point temperature, units = K
    return tdps

def _calc_relhumid(tas, tdps):
    """
    Calculate relative humidity
    Inputs: air temperature, dewpoint temperature (both K)
    Returns: realtive humidity (%, or 0-100 scale)
    """
    es = 0.611 * np.exp(5423 * ((1/273) - (1/tas)))   # calculates saturation vapor pressure using air temp
    e_vap = 0.611 * np.exp(5423 * ((1/273) - (1/tdps)))   # calculates vapor pressure using dew point temp
    hurs = 100 * (e_vap/es)
    return hurs

def _calc_windmag(u10, v10):
    """
    Calculates wind speed
    Inputs: u and v wind components (both m/s)
    Returns: wind speed/magnitude (m/s)
    """
    sfcWind = np.sqrt((u10)**2  + (v10)**2)   # calculates wind magnitude, units = ms-1
    return sfcWind

def _calc_winddir(u10, v10):    # This function would only be needed if wind components are provided, but direction is not
    """
    Calculates wind direction
    Inputs: u and v wind components (both m/s)
    Returns: wind direction (degrees)
    NOTE: this is not functional at the moment due to performance issues
    """
    pass        # this is a complicated calculation -- looking for options
    #return sfcWind_dir

def _calc_ps(psl, elev, temp):
    """
    Calculates station air pressure from sea level air pressure, if station pressure is not available
    Inputs: sea level pressure (Pa), elevation (m), and air temperature (K)
    Returns: air pressure (Pa)
    Note: this calculation checks with this formula, with differences due to rounding in the decimal place:
    https://keisan.casio.com/exec/system/1224575267
    """
    ps = psl / ((1 - ((0.0065 * elev)/(temp + 0.0065 * elev)))**-5.257)
    return ps

def _calc_ps_alt(alt, elev):
    """
    Calculates station air pressure from altimeter setting and station elevation, if station pressure is not available
    Inputs: altimeter setting (Pa) and station elevation (m)
    Returns: air pressure (Pa)
    Note: this calculation uses the following formula:
    https://www.weather.gov/media/epz/wxcalc/stationPressure.pdf
    """
    alt = alt / 3386.39 # Convert altimeter from Pa to inHg for use in formula
    ps = alt * ((288-0.0065*elev)/288)**5.2561
    ps = _unit_pres_inHg_to_pa(ps) # Convert back to Pa from inHg
    return ps
