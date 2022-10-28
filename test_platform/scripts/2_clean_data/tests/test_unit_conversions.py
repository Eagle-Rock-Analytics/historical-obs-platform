"""
This script runs tests that unit conversions perform as expected
"""
import pandas as pd
import pytest
from calc_clean import (_unit_degC_to_K,
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
                        _lon_DMm_to_Dd)

test_data = pd.read_csv("../tests/test_dataset.csv")

## -----------------------------------------------------------------------------------------
## Temperature conversions -- desired unit is K
@pytest.fixture
def grab_temp_C(test_data):
    return test_data['temp_c']

def grab_temp_f(test_data):
    return test_data['temp_f']

def test_temp_conversion_degC(grab_temp_C):
    """Test that the _unit_degC_to_K function correctly converts from degC to K"""
    calc_clean_converted = _unit_degC_to_K(grab_temp_C)
    correct_conversion = (grab_temp_C + 273.15)
    assert correct_conversion.equals(calc_clean_converted)

def test_temp_conversion_degF(grab_temp_f):
    """Test that the _unit_degF_to_K function correctly converts from degF to K"""
    calc_clean_converted = _unit_degF_to_K(grab_temp_f)
    correct_conversion = (5/9) * (grab_temp_f - 32) + 273.15
    assert correct_conversion.equals(calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Precipitation conversions -- desired unit is mm / TIME
## Note: Scaling to hourly data will occur in qa/qc step
## Note: Therefore this unit test for the cleaning stage will ensure that data is in mm per unit of time at present pre scaling
@pytest.fixture
def grab_precip(test_data):
    """Grab the precipitation variable. Native units could be mm or inches"""
    return test_data['precip_in']

def test_precip_conversion_inches(grab_precip):
    """Test that the _unit_precip_in_to_mm function correctly converts from inches to mm"""
    calc_clean_converted = _unit_precip_in_to_mm(grab_precip)
    correct_conversion = grab_precip * 25.4
    assert correct_conversion.equals(calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Wind conversions -- desired unit is m/s
@pytest.fixture
def grab_wind_kts(test_data):
    """Grab the wind variable. Native units could be m/s, mph, knots"""
    return test_data['wind_kts']

def grab_wind_mph(test_data):
    return test_data['wind_mph']

def test_wind_conversion_kts(grab_wind_kts):
    """Test that the _unit_windspd_kts_to_ms correctly converts from knots to m/s"""
    calc_clean_converted = _unit_windspd_kts_to_ms(grab_wind_kts)
    correct_conversion = grab_wind_kts / 1.94
    assert correct_conversion.equals(calc_clean_converted)

def test_wind_conversion_mph(grab_wind_mph):
    """Test that the _unit_windspd_mph_to_ms correctly converts from mph to m/s"""
    calc_clean_converted = _unit_windspd_mph_to_ms(grab_wind_mph)
    correct_conversion = grab_wind_mph / 2.237
    assert correct_conversion.equals(calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Pressure conversions -- desired unit is Pascal
@pytest.fixture
def grab_pressure_hpa(test_data):
    """Grab the pressure variable. Native units could be hPa, mb, Pa, inHg"""
    return test_data['press_hpa']

def grab_pressure_inhg(test_data):
    return test_data['press_inhg']

def test_pressure_conversion_hpa(grab_pressure_hpa):
    """Test that the _unit_pres_hpa_to_pa correctly converts from hPa to Pa"""
    calc_clean_converted = _unit_pres_hpa_to_pa(grab_pressure_hpa)
    correct_conversion = grab_pressure_hpa * 100.
    assert correct_conversion.equals(calc_clean_converted)

def test_pressure_conversion_inHg(grab_pressure_inhg):
    """Test that the _unit_pres_inHg_to_pa correctly converts from inHg to Pa"""
    calc_clean_converted = _unit_pres_inHg_to_pa(grab_pressure_inhg)
    correct_conversion = grab_pressure_inhg * 3386.39
    assert correct_conversion.equals(calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Moisture conversions -- desired unit is kg/kg
## Note: Not a primary variable, but optionally needed for some derivations
@pytest.fixture
def grab_moisture(test_data):
    """Grabs the water vapor mixing ratio. Native units could be g/kg or kg/kg"""
    return test_data['moisture_kgkg']

def test_moisture_conversion_kgkg(grab_moisture):
    """Test that the _unit_moisture_gkg_to_kgkg correctly converts from g/kg to kg/kg"""
    calc_clean_converted = _unit_moisture_gkg_to_kgkg(grab_moisture)
    correct_conversion = grab_moisture / 1000.
    assert correct_conversion.equals(calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Elevation conversions -- desired unit is m
@pytest.fixture
def grab_elevation(test_data):
    """Grabs the elevation. Native units could be m or feet"""
    return test_data['elev_feet'] ## check variable name

def test_elevation_conversion_feet(grab_elevation):
    """Test that the _unit_elev_ft_to_m correctly converts from feet to m"""
    calc_clean_converted = _unit_elev_ft_to_m(grab_elevation)
    correct_conversion = grab_elevation * 0.3048
    assert correct_conversion.equals(calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Latitude conversions -- desired unit is decimal degrees
@pytest.fixture
def grab_latitude_dms(test_data):
    """Grabs the latitude. Native units could be decimal-minutes-seconds (DMS), decimal degrees (D.d), or LORAN (DM.m)"""
    return test_data['lat_dms'] ## check variable name

def grab_latitude_dmm(test_data):
    return test_data['lat_dmm']

def test_latitude_conversion_dms(grab_latitude_dms):
    """Test that the _lat_dms_to_dd correctly converts from DMS to D.d"""
    calc_clean_converted = _lat_dms_to_dd(grab_latitude_dms)
    correct_conversion = float(grab_latitude_dms[:2]) + float(grab_latitude_dms[3:5])/60 + float(grab_latitude_dms[6:])/3600
    assert correct_conversion.equals(calc_clean_converted)

def test_latitude_conversion_dmm(grab_latitude_dmm):
    """Test that the _lat_DMm_to_Dd correctly converts from DMm to D.d"""
    calc_clean_converted = _lat_DMm_to_Dd(grab_latitude_dmm)
    correct_conversion = float(grab_latitude_dmm[:2]) + float(grab_latitude_dmm[2:])/60
    assert correct_conversion.equals(calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Longitude conversions -- desired unit is decimal degrees
@pytest.fixture
def grab_longitude_dms(test_data):
    """Grabs the longitude. Native units could be decimal-minutes-seconds (DMS), decimal degrees (D.d), or LORAN (DM.m)"""
    return test_data['lon_dms'] ## check variable name

def grab_longitude_dmm(test_data):
    return test_data['lon_dmm']

def test_longitude_conversion_dms(grab_longitude_dms):
    """Test that the _lon_dms_to_dd correctly converts from DMS to D.d"""
    calc_clean_converted = _lon_dms_to_dd(grab_longitude_dms)

    if grab_longitude_dms[0] == "-":
        grab_longitude_dms.strip('-')
    correct_conversion = -1 * (float(grab_longitude_dms[:3]) + (float(grab_longitude_dms[4:6])/60) + (float(grab_longitude_dms[7:])/3600))
    assert correct_conversion.equals(calc_clean_converted)

def test_longitude_conversion_dmm(grab_longitude_dmm):
    """Test that the _lon_DMm_to_Dd correctly converts from DM.m to D.d"""
    calc_clean_converted = _lon_DMm_to_Dd(grab_longitude_dmm)
    correct_conversion = -1 * (float(grab_longitude_dmm[:3]) + (float(grab_longitude_dmm[3:])/60))
    assert correct_conversion.equals(calc_clean_converted)
