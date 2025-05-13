"""
This script performs merge protocols for deriving any missing variables for ingestion into the Historical Observations Platform, 
and is independent of network. Missing variables are defined as variables that can be calculated for which there are the 
required sensors to calculate a variable. Observed and calculated data are not mixed, i.e., the missing variable derivation 
is not meant to in-fill potential valid missing observations.

Variables that can be derived if observations are missing
- Dewpoint temperature: requires air temperature, relative humidity; or vapor pressure (not kept)
- Relative humidity: requires air temperature, dewpoint temperature
- Pressure: air temperature, elevation, sea level pressure
- Air temperature: requires dewpoint temperature, relative humidity

Variables that cannot be derived if observations are missing
- Wind speed, direction: requires u and v components (not kept)
- Surface radiation
- Precipitation
"""

# Import libraries
import pandas as pd
import numpy as np
import xarray as xr
import logging

# Import relevant merge util functions
from merge_utils import get_file_paths

# why breaking????


## Identify vars that can be derived
def merge_derive_missing_vars(df: pd.DataFrame) -> pd.DataFrame:
    """
    
    Parameters
    ----------
    df : pd.DataFrame
        input dataframe

    Returns
    -------
    df : pd.DataFrame
        dataframe with newly added derived variable
    """

    # vars that can be derived
    derive_vars = ['tdps', 'hurs', 'ps', 'tas'] # pressure is funky; what about tdps_derived

    vars_to_check = df.columns

    # first check if station has any vars that can be derived
    for item in derive_vars:
        if item in df.columns:
            print(f"{item} is available!")
            continue

        else: # var is missing
            # check if required inputs are available
            if item == "tdps":
                print('hello dewpoint')

    return None

# --------------------------------------------------
def input_var_check(df: pd.DataFrame, var1: str, var2: str, var3: str | None = None) -> bool:
    """
    Flexible check if required secondary input variables are available to derive a primary variable.

    Parameters
    ----------
    df : pd.DataFrame
        input dataframe to check against
    var1 : str
        name of secondary input var 1
    var2 : str
        name of secondary input var 2
    var3 : str
        name of secondary input var 3; optional, only applicable for air pressure

    Returns
    -------
    bool
        True if all required input vars present; False if not
    """

    if var3 == None:
        if var1 in df.columns and var2 in df.columns:
            return True
        else: 
            return False
    else:
        if var1 in df.columns and var2 in df.columns and var3 in df.columns:
            return True
        else:
            return False


## Derived variable calculations
def _calc_dewpointtemp(tas: pd.Series, hurs: pd.Series) -> pd.Series:
    """Calculates dew point temperature, method 1

    Parameters
    ----------
    tas : pd.Series
        air temperature, K
    hurs: pd.Series
        relative humidity, % or 0-100

    Returns
    -------
    tdps : pd.Series
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


def _calc_relhumid(tas: pd.Series, tdps: pd.Series) -> pd.Series:
    """Calculate relative humidity

    Parameters
    ----------
    tas : pd.Series
        air temperature, K
    tdps : pd.Series
        dewpoint temperature, K

    Returns
    -------
    hurs : pd.Series
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


def _calc_ps(psl: pd.Series, elev: pd.Series, temp: pd.Series) -> pd.Series:
    """Calculates station air pressure from sea level air pressure, if station pressure is not available

    Parameters
    -----------
    psl : pd.Series
        sea level pressure, Pa
    elev : pd.Series
        elevation, m
    temp : pd.Series
        air temperature, K

    Returns
    -------
    ps : pd.Series
        surface pressure, Pa

    Notes
    -----
    1. This calculation checks with this formula, with differences due to rounding in the decimal place:
    https://keisan.casio.com/exec/system/1224575267
    """
    ps = psl / ((1 - ((0.0065 * elev) / (temp + 0.0065 * elev))) ** -5.257)
    return ps


def _calc_ps_alt(alt: pd.Series, elev: pd.Series) -> pd.Series:
    """Calculates station air pressure from altimeter setting and station elevation, if station pressure is not available

    Parameters
    ----------
    alt : pd.Series
        altimeter setting, Pa
    elev : pd.Series
        elevation, m

    Returns
    -------
    ps : pd.Series
        air pressure, Pa

    References
    ----------
    [1] This calculation uses the following formula: https://www.weather.gov/media/epz/wxcalc/stationPressure.pdf
    """

    def _unit_pres_inHg_to_pa(data: pd.Series) -> pd.Series:
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

    alt = alt / 3386.39  # Convert altimeter from Pa to inHg for use in formula
    ps = alt * ((288 - 0.0065 * elev) / 288) ** 5.2561
    ps = _unit_pres_inHg_to_pa(ps)  # Convert back to Pa from inHg
    return ps