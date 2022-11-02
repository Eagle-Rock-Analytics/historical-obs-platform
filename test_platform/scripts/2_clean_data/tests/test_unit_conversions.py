"""
This script runs tests that unit conversions perform as expected.

To run: navigate to the folder where this is stored, and run "pytest test_unit_conversions.py" or just "pytest" from the command line
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
                        _unit_moisture_gkg_to_kgkg,
                        _lat_dms_to_dd,
                        _lat_DMm_to_Dd,
                        _lon_dms_to_dd,
                        _lon_DMm_to_Dd)

@pytest.fixture
def df():
    return pd.read_csv('test_dataset.csv')

## -----------------------------------------------------------------------------------------
## Temperature conversions -- desired unit is K
@pytest.fixture
def grab_temp_C(df):
    return df['temp_c']

@pytest.fixture
def grab_temp_f(df):
    return df['temp_f']

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
def grab_precip(df):
    """Grab the precipitation variable. Native units could be mm or inches"""
    return df['precip_in']

def test_precip_conversion_inches(grab_precip):
    """Test that the _unit_precip_in_to_mm function correctly converts from inches to mm"""
    calc_clean_converted = _unit_precip_in_to_mm(grab_precip)
    correct_conversion = grab_precip * 25.4
    assert correct_conversion.equals(calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Wind conversions -- desired unit is m/s
@pytest.fixture
def grab_wind_kts(df):
    """Grab the wind variable. Native units could be m/s, mph, knots"""
    return df['wind_kts']

@pytest.fixture
def grab_wind_mph(df):
    return df['wind_mph']

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
def grab_pressure_hpa(df):
    """Grab the pressure variable. Native units could be hPa, mb, Pa, inHg"""
    return df['press_hpa']

@pytest.fixture
def grab_pressure_inhg(df):
    return df['press_inhg']

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
def grab_moisture(df):
    """Grabs the water vapor mixing ratio. Native units could be g/kg or kg/kg"""
    return df['moisture_gkg']

def test_moisture_conversion_kgkg(grab_moisture):
    """Test that the _unit_moisture_gkg_to_kgkg correctly converts from g/kg to kg/kg"""
    calc_clean_converted = _unit_moisture_gkg_to_kgkg(grab_moisture)
    correct_conversion = grab_moisture / 1000.
    assert correct_conversion.equals(calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Elevation conversions -- desired unit is m
@pytest.fixture
def grab_elevation(df):
    """Grabs the elevation. Native units could be m or feet"""
    return df['elev_feet']

def test_elevation_conversion_feet(grab_elevation):
    """Test that the _unit_elev_ft_to_m correctly converts from feet to m"""
    calc_clean_converted = _unit_elev_ft_to_m(grab_elevation)
    correct_conversion = grab_elevation * 0.3048
    assert correct_conversion.equals(calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Latitude conversions -- desired unit is decimal degrees
@pytest.fixture
def grab_latitude_dms(df):
    """Grabs the latitude. Native units could be decimal-minutes-seconds (DMS), decimal degrees (D.d), or LORAN (DM.m)"""
    return df['lat_dms'] ## check variable name

@pytest.fixture
def grab_latitude_dmm(df):
    return df['lat_dmm']

def test_latitude_conversion_dms(grab_latitude_dms):
    """Test that the _lat_dms_to_dd correctly converts from DMS to D.d"""
    for i in range(len(grab_latitude_dms)):
        calc_clean_converted = _lat_dms_to_dd(grab_latitude_dms[i])
        _deg = float(grab_latitude_dms[i][:2])
        _min = float(grab_latitude_dms[i][3:5])
        _sec = float(grab_latitude_dms[i][6:])
    correct_conversion = _deg + _min/60 + _sec/3600
    assert correct_conversion == calc_clean_converted

# def test_latitude_conversion_dmm(grab_latitude_dmm):
#     """Test that the _lat_DMm_to_Dd correctly converts from DMm to D.d"""
#     for i in range(len(grab_latitude_dmm)):
#         calc_clean_converted = _lat_DMm_to_Dd(grab_latitude_dmm[i])
#         _lat = grab_latitude_dmm['lat_dmm'][:2].apply(lambda x: float(x))
#         _min = grab_latitude_dmm['lat_dmm'][2:].apply(lambda y: float(y))
#     correct_conversion = _lat + _min/60
#     assert correct_conversion == (calc_clean_converted)


## -----------------------------------------------------------------------------------------
## Longitude conversions -- desired unit is decimal degrees
@pytest.fixture
def grab_longitude_dms(df):
    """Grabs the longitude. Native units could be decimal-minutes-seconds (DMS), decimal degrees (D.d), or LORAN (DM.m)"""
    return df['lon_dms'] ## check variable name

@pytest.fixture
def grab_longitude_dmm(df):
    return df['lon_dmm']

def test_longitude_conversion_dms(grab_longitude_dms):
    """Test that the _lon_dms_to_dd correctly converts from DMS to D.d"""
    for i in range(len(grab_longitude_dms)):
        calc_clean_converted = _lon_dms_to_dd(grab_longitude_dms[i])
        if grab_longitude_dms[i][0] == "-":
            grab_longitude_dms.strip('-')

    correct_conversion = -1 * (float(grab_longitude_dms[i][:3]) + (float(grab_longitude_dms[i][4:6])/60) + (float(grab_longitude_dms[i][7:])/3600))
    assert correct_conversion == (calc_clean_converted)

# def test_longitude_conversion_dmm(grab_longitude_dmm):
#     """Test that the _lon_DMm_to_Dd correctly converts from DM.m to D.d"""
#     for i in range(len(grab_longitude_dmm)):
#         calc_clean_converted = _lon_DMm_to_Dd(grab_longitude_dmm[i])
#
#     correct_conversion = -1 * (float(grab_longitude_dmm[i][:3]) + (float(grab_longitude_dmm[i][3:])/60))
#     assert correct_conversion == (calc_clean_converted)
