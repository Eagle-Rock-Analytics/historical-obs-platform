"""
This is a script where Stage 2: Clean related common functions, conversions, and operations is stored for ease of use
for the Historical Observations Platform.
"""

## Process to call this script: from SCRIPT import FUNCTION
## Example: from calc_pull import get_wecc_poly

## Import Libraries
import geopandas as gp
import numpy as np

##---------------------------------------------------------------------------------------------
## Useful functions
def get_wecc_poly(terrpath, marpath):
    """Identifies a bbox of WECC area to filter stations against.

    Parameters
    ----------
    terrpath : str
        shapefile for terrestrial WECC boundary
    marpath : str
        shapefile for maritime WECC boundary
    
    Returns
    -------
    t : gp.Polygon
        polygon of terrestrial WECC boundary
    m : gp.Polygon
        polygon of maritime WECC boundary
    bbox : gp.Polygon
        bounding box of t and m union
    """
    t = gp.read_file(terrpath)  ## Read in terrestrial WECC shapefile.
    m = gp.read_file(marpath)  ## Read in marine WECC shapefile.
    bbox = t.union(m).bounds  ## Combine polygons and get bounding box of union.
    return t, m, bbox

##---------------------------------------------------------------------------------------------
## Unit conversions, if required
def _unit_degC_to_K(data):
    """Converts temperature from degC to K

    Parameters
    ----------
    data : float
        input data to convert
    
    Returns
    -------
    data : float
        data converted to K
    """
    data = data + 273.15
    return data

##---------------------------------------------------------------------------------------------
def _unit_degF_to_K(data):
    """Converts temperature from degF to K
    Parameters
    ----------
    data : float
        input data to convert
    
    Returns
    -------
    data : float
        data converted to K
    """
    data = (5 / 9) * (data - 32) + 273.15
    return data

##---------------------------------------------------------------------------------------------
def _unit_pres_hpa_to_pa(data):
    """Converts air pressure from hectopascals to pascals

    Parameters
    ----------
    data : float
        input data to convert
    
    Returns
    -------
    data : float
        data converted to Pa

    Notes
    ------
    This also works for the conversion from mb
    """
    data = data * 100.0
    return data

##---------------------------------------------------------------------------------------------
def _unit_pres_kpa_to_pa(data):
    """Converts air pressure from kilopascals to pascals

    Parameters
    ----------
    data : float
        input data to convert
    
    Returns
    -------
    data : float
        data converted to Pa
    """
    data = data * 1000.0
    return data

##---------------------------------------------------------------------------------------------
def _unit_pres_inHg_to_pa(data):
    """Converts air pressure from inHg to hectopascals
    
    Parameters
    ----------
    data : float
        input data to convert
    
    Returns
    -------
    data : float
        data converted to Pa
    """
    data = data * 3386.39
    return data

##---------------------------------------------------------------------------------------------
## Wind speed conversions: Desired working unit should be m/s
def _unit_windspd_kts_to_ms(data):
    """Converts windspeed from knots to m/s

    Parameters
    ----------
    data : float
        input data to convert
    
    Returns
    -------
    data : float
        data converted to m/s
    """
    data = data / 1.94
    return data

##---------------------------------------------------------------------------------------------
def _unit_windspd_mph_to_ms(data):
    """Converts windspeed from miles-per-hour to m/s

    Parameters
    ----------
    data : float
        input data to convert
    
    Returns
    -------
    data : float
        data converted to m/s
    """
    data = data / 2.237
    return data

##---------------------------------------------------------------------------------------------
## Moisture conversions: Desired working unit should be  kg/kg
def _unit_moisture_gkg_to_kgkg(data):
    """Converts moisture ratios from g/kg to kg/kg

    Parameters
    ----------
    data : float
        input data to convert
    
    Returns
    -------
    data : float
        data converted to kg/kg
    """
    data = data / 1000.0
    return data

##---------------------------------------------------------------------------------------------
## Precipitation conversions: Desired working unit should be mm/TIME
## Note on desired unit - the conversion to hourly as a rate will occur in calc_qaqc
def _unit_precip_in_to_mm(data):
    """Converts precipitation from inches to mm

    Parameters
    ----------
    data : float
        input data to convert
    
    Returns
    -------
    data : float
        data converted to mm
    """
    data = data * 25.4
    return data

##---------------------------------------------------------------------------------------------
## Elevation conversions: Desired working unit should be meters
def _unit_elev_ft_to_m(data):
    """Converts elevation from feet to meters

    Parameters
    ----------
    data : float
        input data to convert
    
    Returns
    -------
    data : float
        data converted to m
    """
    data = data * 0.3048
    return data

##---------------------------------------------------------------------------------------------
## Latitude/Longitude conversions: Desired working unit should be decimal degrees N/W
def _lat_dms_to_dd(data):
    """Converts latitude from decimal-minutes-seconds to decimal degrees

    Parameters
    ----------
    data : float
        input data to convert
    
    Returns
    -------
    data : float
        data converted to dd
    """
    data = float(data[:2]) + float(data[3:5]) / 60 + float(data[6:]) / 3600
    return data

##---------------------------------------------------------------------------------------------
def _lon_dms_to_dd(data):
    """Converts longitude from decimal-minutes-seconds to decimal degrees 
    and ensures that western hemisphere lons are negative by convention

    Parameters
    ----------
    data : float
        input data to convert
    
    Returns
    -------
    data : float
        data converted to dd
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

##---------------------------------------------------------------------------------------------
def _lon_DMm_to_Dd(data):
    """This is specific to CWOP longitude data converting from LORAN (DM.m) coordinates to decimal-degrees (D.d) for the WESTERN HEMISPHERE.

    Parameters
    ----------
    data : float
        input data to convert
    
    Returns
    -------
    data : float
        data converted to D.d
    """
    _min = float(data[:3])
    _sec = float(data[3:])
    data = -1 * (_min + _sec / 60)
    return data

##---------------------------------------------------------------------------------------------
def _lat_DMm_to_Dd(data):
    """This is specific to CWOP latitude data converting from LORAN (DM.m) coordinates to decimal-degrees (D.d).

    Parameters
    ----------
    data : float
        input data to convert
    
    Returns
    -------
    data : float
        data converted to D.d
    """
    _deg = float(data[:2])
    _mm = float(data[2:])
    data = _deg + _mm / 60
    return data


##---------------------------------------------------------------------------------------------
## Derived variable calculations
def _calc_dewpointtemp_opt1(tas, hurs):
    """Calculates dew point temperature, method 1

    Parameters
    ----------
    tas : float
        air temperature, K
    hurs: float
        relative humidity, % or 0-100
    
    Returns
    -------
    tdps : float
        dewpoint temperature, K
    """
    es = 0.611 * np.exp(
        5423 * ((1 / 273) - (1 / tas))
    )  # calculates saturation vapor pressure
    e_vap = (
        es * hurs
    ) / 100.0  # calculates vapor pressure, IF NOT ALREADY OBSERVED -- will need ifelse statement
    tdps = (
        (1 / 273) - 0.0001844 * np.log(e_vap / 0.611)
    ) ** -1  # calculates dew point temperature, units = K
    return tdps

##---------------------------------------------------------------------------------------------
def _calc_dewpointtemp_opt2(e_vap):
    """Calculates dew point temperature, method 2

    Parameters
    ----------
    e_vap : float
        vapor pressure, Pa
    
    Returns
    -------
    tdps : float
        dewpoint temperature, K
    """
    tdps = (
        (1 / 273) - 0.0001844 * np.log(e_vap / 0.611)
    ) ** -1  # calculates dew point temperature, units = K
    return tdps

##---------------------------------------------------------------------------------------------
def _calc_relhumid(tas, tdps):
    """Calculate relative humidity

    Parameters
    ----------
    tas : float
        air temperature, K
    tdps : float
        dewpoint temperature, K
    
    Returns
    -------
    hurs : float
        relative humidity, % (0-100)
    """

    es = 0.611 * np.exp(
        5423 * ((1 / 273) - (1 / tas))
    )  # calculates saturation vapor pressure using air temp
    e_vap = 0.611 * np.exp(
        5423 * ((1 / 273) - (1 / tdps))
    )  # calculates vapor pressure using dew point temp
    hurs = 100 * (e_vap / es)
    return hurs

##---------------------------------------------------------------------------------------------
def _calc_windmag(u10, v10):
    """Calculates wind speed

    Parameters
    ----------
    u10 : float
        u-direction wind
    v10 : float
        v-direction wind   
    
    Returns
    -------
    sfcWind : float
        wind magnitude, m/s

    Notes
    -----
    1. u and v wind components (both in m/s)
    """
    sfcWind = np.sqrt(
        (u10) ** 2 + (v10) ** 2
    )  # calculates wind magnitude, units = ms-1
    return sfcWind

##---------------------------------------------------------------------------------------------
def _calc_winddir(
    u10, v10
):  
    """Calculates wind direction

    Parameters
    ----------
    u10 : float
        u-direction wind
    v10 : float
        v-direction wind   

    Returns
    -------
    sfcWind_dir : float
        wind direction, degrees

    Notes
    -----
    1. u and v wind components (both in m/s)
    2. This is NOT functional at the moment due to performance issues
    3. This function would only be needed if wind components are provided, but direction is not
    """
    pass  # this is a complicated calculation -- looking for options
    # return sfcWind_dir

##---------------------------------------------------------------------------------------------
def _calc_ps(psl, elev, temp):
    """Calculates station air pressure from sea level air pressure, if station pressure is not available

    Parameters
    -----------
    psl : float
        sea level pressure, Pa
    elev : float
        elevation, m
    temp : float
        air temperature, K

    Returns
    -------
    ps : float
        surface pressure, Pa

    Notes
    ----- 
    1. This calculation checks with this formula, with differences due to rounding in the decimal place:
    https://keisan.casio.com/exec/system/1224575267
    """
    ps = psl / ((1 - ((0.0065 * elev) / (temp + 0.0065 * elev))) ** -5.257)
    return ps

##---------------------------------------------------------------------------------------------
def _calc_ps_alt(alt, elev):
    """Calculates station air pressure from altimeter setting and station elevation, if station pressure is not available

    Parameters
    ----------
    alt : float
        altimeter setting, PA
    elev : float
        elevation, m
    
    Returns
    -------
    ps : float
        air pressure, Pa

    Notes
    -----
    This calculation uses the following formula: https://www.weather.gov/media/epz/wxcalc/stationPressure.pdf
    """
    alt = alt / 3386.39  # Convert altimeter from Pa to inHg for use in formula
    ps = alt * ((288 - 0.0065 * elev) / 288) ** 5.2561
    ps = _unit_pres_inHg_to_pa(ps)  # Convert back to Pa from inHg
    return ps
