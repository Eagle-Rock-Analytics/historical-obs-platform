"""
This script runs tests that unit conversions perform as expected
"""

import pytest
from calc_clean import _unit_degC_to_K,
                       _unit_degF_to_K,
                       _unit_precip_in_to_mm,
                       _unit_windspd_kts_to_ms,
                       _unit_windspd_mph_to_ms,
                       _unit_pres_hpa_to_pa,
                       _unit_pres_inHg_to_pa,
                       _unit_elev_ft_to_m,
                       _lat_dms_to_dd,
                       _lat_DMm_to_Dd,
                       _lon_dms_to_dd,
                       _lon_DMm_to_Dd


## -----------------------------------------------------------------------------------------
## Temperature conversions -- desired unit is K
@pytest.fixture
def grab_temp():
    """Grab the temperature variable. Native unit could be K, degC, or degF"""
    return temp
    # Note, dewpoint temperature is also appropriately used here

def test_temp_conversion_degC(grab_temp):
    """Test that the _unit_degC_to_K function correctly converts from degC to K"""
    calc_clean_converted = _unit_degC_to_K(grab_temp)
    correct_conversion = (grab_temp + 273.15)
    assert correct_conversion.equals(calc_clean_converted)

def test_temp_conversion_degF(grab_temp):
    """Test that the _unit_degF_to_K function correctly converts from degF to K"""
    calc_clean_converted = _unit_degF_to_K(grab_temp)
    correct_conversion = (5/9) * (grab_temp - 32) + 273.15
    assert correct_conversion.equals(calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Precipitation conversions -- desired unit is mm / TIME
## Note: Scaling to hourly data will occur in qa/qc step
## Note: Therefore this unit test for the cleaning stage will ensure that data is in mm per unit of time at present pre scaling
@pytest.fixture
def grab_precip():
    """Grab the precipitation variable. Native units could be mm or inches"""
    return pr

def test_precip_conversion_inches(grab_precip):
    """Test that the _unit_precip_in_to_mm function correctly converts from inches to mm"""
    calc_clean_converted = _unit_precip_in_to_mm(grab_precip)
    correct_conversion = grab_precip * 25.4
    assert correct_conversion.equals(calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Wind conversions -- desired unit is m/s
@pytest.fixture
def grab_wind():
    """Grab the wind variable. Native units could be m/s, mph, knots"""
    return sfcWind

def test_wind_conversion_kts(grab_wind):
    """Test that the _unit_windspd_kts_to_ms correctly converts from knots to m/s"""
    calc_clean_converted = _unit_windspd_kts_to_ms(grab_wind)
    correct_conversion = grab_wind / 1.94
    assert correct_conversion.equals(calc_clean_converted)

def test_wind_conversion_mph(grab_wind):
    """Test that the _unit_windspd_mph_to_ms correctly converts from mph to m/s"""
    calc_clean_converted = _unit_windspd_mph_to_ms(grab_wind)
    correct_conversion = grab_wind / 2.237
    assert correct_conversion.equals(calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Pressure conversions -- desired unit is Pascal
@pytest.fixture
def grab_pressure():
    """Grab the pressure variable. Native units could be hPa, mb, Pa, inHg"""
    return ps

def test_pressure_conversion_hpa(grab_pressure):
    """Test that the _unit_pres_hpa_to_pa correctly converts from hPa to Pa"""
    calc_clean_converted = _unit_pres_hpa_to_pa(grab_pressure)
    correct_conversion = grab_pressure * 100.
    assert correct_conversion.equals(calc_clean_converted)

def test_pressure_conversion_inHg(grab_pressure):
    """Test that the _unit_pres_inHg_to_pa correctly converts from inHg to Pa"""
    calc_clean_converted = _unit_pres_inHg_to_pa(grab_pressure)
    correct_conversion = grab_pressure * 3386.39
    assert correct_conversion.equals(calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Moisture conversions -- desired unit is kg/kg
## Note: Not a primary variable, but optionally needed for some derivations
@pytest.fixture
def grab_moisture():
    """Grabs the water vapor mixing ratio. Native units could be g/kg or kg/kg"""
    return r ## check this variable name

def test_moisture_conversion_kgkg(grab_moisture):
    """Test that the _unit_moisture_gkg_to_kgkg correctly converts from g/kg to kg/kg"""
    calc_clean_converted = _unit_moisture_gkg_to_kgkg(grab_moisture)
    correct_conversion = grab_moisture / 1000.
    assert correct_conversion.equals(calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Elevation conversions -- desired unit is m
@pytest.fixture
def grab_elevation():
    """Grabs the elevation. Native units could be m or feet"""
    return elev ## check variable name

def test_elevation_conversion_feet(grab_elevation):
    """Test that the _unit_elev_ft_to_m correctly converts from feet to m"""
    calc_clean_converted = _unit_elev_ft_to_m(grab_elevation)
    correct_conversion = grab_elevation * 0.3048
    assert correct_conversion.equals(calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Latitude conversions -- desired unit is decimal degrees
@pytest.fixture
def grab_latitude():
    """Grabs the latitude. Native units could be decimal-minutes-seconds (DMS), decimal degrees (D.d), or LORAN (DM.m)"""
    return lat ## check variable name

def test_latitude_conversion_dms(grab_latitude):
    """Test that the _lat_dms_to_dd correctly converts from DMS to D.d"""
    calc_clean_converted = _lat_dms_to_dd(grab_latitude)
    correct_conversion = float(grab_latitude[:2]) + float(grab_latitude[3:5])/60 + float(grab_latitude[6:])/3600
    assert correct_conversion.equals(calc_clean_converted)

def test_latitude_conversion_dmm(grab_latitude):
    """Test that the _lat_DMm_to_Dd correctly converts from DMm to D.d"""
    calc_clean_converted = _lat_DMm_to_Dd(grab_latitude)
    correct_conversion = float(grab_latitude[:2]) + float(grab_latitude[2:])/60
    assert correct_conversion.equals(calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Longitude conversions -- desired unit is decimal degrees
@pytest.fixture
def grab_longitude():
    """Grabs the longitude. Native units could be decimal-minutes-seconds (DMS), decimal degrees (D.d), or LORAN (DM.m)"""
    return lon ## check variable name

def test_longitude_conversion_dms(grab_longitude):
    """Test that the _lon_dms_to_dd correctly converts from DMS to D.d"""
    calc_clean_converted = _lon_dms_to_dd(grab_longitude)

    if grab_longitude[0] == "-":
        grab_longitude.strip('-')
    correct_conversion = -1 * (float(grab_longitude[:3]) + (float(grab_longitude[4:6])/60) + (float(grab_longitude[7:])/3600))
    assert correct_conversion.equals(calc_clean_converted)

def test_longitude_conversion_dmm(grab_longitude):
    """Test that the _lon_DMm_to_Dd correctly converts from DM.m to D.d"""
    calc_clean_converted = _lon_DMm_to_Dd(grab_longitude)
    correct_conversion = -1 * (float(grab_longitude[:3]) + (float(grab_longitude[3:])/60))
    assert correct_conversion.equals(calc_clean_converted)
