"""
This script runs tests that variable derivations perform as expected.

To run: type "pytest test_derived_variables.py" from the command line
"""
import pandas as pd
import pytest
import sys
import numpy as np

## Adding 2_clean_data.calc_clean to system path
## Note to use: must insert full path to where the 2_clean_data dir is stored locally
sys.path.insert(0, '/path/to/local/historical-obs/historical-obs-platform/test_platform/scripts/2_clean_data/')

from calc_clean import (_calc_dewpointtemp_opt1,
                        _calc_dewpointtemp_opt2,
                        _calc_relhumid,
                        _calc_windmag,
                        _calc_winddir,
                        _calc_ps,
                        _calc_ps_alt,
                        _unit_pres_inHg_to_pa)

@pytest.fixture
def df():
    return pd.read_csv('test_dataset.csv')

## -----------------------------------------------------------------------------------------
## Dewpoint temperature
@pytest.fixture
def grab_temp_k(df):
    return df['temp_k']

@pytest.fixture
def grab_hurs(df):
    return df['hurs']

@pytest.fixture
def grab_vapor_pres_pa(df):
    return df['vapor_pres_pa']

def test_dewpoint_opt1_derivation(grab_temp_k, grab_hurs):
    """Test that the _calc_dewpointtemp_opt1 function correctly calculates dew point temperature"""
    calc_clean_derived = _calc_dewpointtemp_opt1(grab_temp_k, grab_hurs)
    es = 0.611 * np.exp(5423 * ((1/273) - (1/grab_temp_k)))
    e_vap = (es * grab_hurs) / 100.
    correct_derivation = ((1/273) - 0.0001844 * np.log(e_vap/0.611))**-1
    assert correct_derivation.equals(calc_clean_derived)

def test_dewpoint_opt2_derivation(grab_vapor_pres_pa):
    """Test that the _calc_dewpointtemp_opt2 function correctly calculates dew point temperature"""
    calc_clean_derived = _calc_dewpointtemp_opt2(grab_vapor_pres_pa)
    correct_derivation = ((1/273) - 0.0001844 * np.log(grab_vapor_pres_pa/0.611))**-1
    assert correct_derivation.equals(calc_clean_derived)

## -----------------------------------------------------------------------------------------
## Relative humidity
@pytest.fixture
def grab_temp_k(df):
    return df['temp_k']

@pytest.fixture
def grab_dewpttemp_k(df):
    return df['dewpoint_temp_k']

def test_relhumid_derivation(grab_temp_k, grab_dewpttemp_k):
    """Test that the _calc_relhumid function correctly calculates relative humidity"""
    calc_clean_derived = _calc_relhumid(grab_temp_k, grab_dewpttemp_k)
    es = 0.611 * np.exp(5423 * ((1/273) - (1/grab_temp_k)))
    e_vap = 0.611 * np.exp(5423 * ((1/273) - (1/grab_dewpttemp_k)))
    correct_derivation = 100. * (e_vap / es)
    assert correct_derivation.equals(calc_clean_derived)

## -----------------------------------------------------------------------------------------
## Wind magnitude
@pytest.fixture
def grab_uwind(df):
    return df['u10']

@pytest.fixture
def grab_vwind(df):
    return df['v10']

def test_windmag_derivation(grab_uwind, grab_vwind):
    """Test that the _calc_windmag function correctly calculates wind magnitude/speed"""
    calc_clean_derived = _calc_windmag(grab_uwind, grab_vwind)
    correct_derivation = np.sqrt((grab_uwind)**2 + (grab_vwind)**2)
    assert correct_derivation.equals(calc_clean_derived)

## -----------------------------------------------------------------------------------------
## Wind direction -- note this is not functional at present

## -----------------------------------------------------------------------------------------
## Station air pressure
@pytest.fixture
def grab_elev(df):
    return df['elev_m']

@pytest.fixture
def grab_temp_k(df):
    return df['temp_k']

@pytest.fixture
def grab_psl_pa(df):
    return df['psl_pa']

@pytest.fixture
def grab_alt(df):
    return df['alt_pa']

def test_pres_derivation_slp(grab_psl_pa, grab_elev, grab_temp_k):
    """Test that the _calc_ps function correctly calculates station air pressure"""
    calc_clean_derived = _calc_ps(grab_psl_pa, grab_elev, grab_temp_k)
    correct_derivation = grab_psl_pa / ((1 - ((0.0065 * grab_elev)/(grab_temp_k + 0.0065 * grab_elev)))**-5.257)
    assert correct_derivation.equals(calc_clean_derived)

def test_pres_derivation_alt(grab_alt, grab_elev):
    """Test that the _calc_ps_alt function correctly calculates station air pressure"""
    calc_clean_derived = _calc_ps_alt(grab_alt, grab_elev)
    alt_inhg = grab_alt / 3386.39
    ps = alt_inhg * ((288-0.0065 * grab_elev)/288)**5.2561
    correct_derivation = _unit_pres_inHg_to_pa(ps)
    assert correct_derivation.equals(calc_clean_derived)
