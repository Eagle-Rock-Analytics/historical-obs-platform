"""
This script performs merge protocols for deriving any missing variables for ingestion into the Historical Observations Platform, 
and is independent of network. Missing variables are defined as variables that can be calculated for which there are the 
required sensors to calculate a variable. Observed and calculated data are not mixed, i.e., the missing variable derivation 
is not meant to in-fill potential valid missing observations.

Variables that can be derived if observations are missing
- Dewpoint temperature: requires air temperature, relative humidity; or vapor pressure (not kept)
- Relative humidity: requires air temperature, dewpoint temperature
- Air temperature: requires dewpoint temperature, relative humidity

Variables that cannot be derived if observations are missing
- Wind speed, direction: requires u and v components (not kept)
- Surface radiation
- Precipitation
- Pressure: air temperature, elevation, sea level pressure (not ideal, tricky calculation)
"""

# Import libraries
import pandas as pd
import numpy as np
import xarray as xr
import logging


## Identify vars that can be derived
def merge_derive_missing_vars(df: pd.DataFrame) -> pd.DataFrame:
    """
    Identifies if any variables can be derived with other input variables.
    If success, variable is derived in the correct unit, attribtues are updated,
    and any flags from the input variables are synergistically flagged.
    If failure, variable is not derived.

    Parameters
    ----------
    df : pd.DataFrame
        input dataframe

    Returns
    -------
    df_flag : pd.DataFrame
        dataframe with newly added derived variable
    """
    print("Running merge_derive_missing_vars...")  # conver to logger when ready

    # vars that can be derived
    derive_vars = ["tdps", "hurs", "tas"]  # only tdps, not tdps_derived

    # first check if station has any vars that can be derived
    for item in derive_vars:
        if item in df.columns:
            print(
                f"{item} is present in station, no derivation necessary."
            )  # convert to logger when set-up
            continue

        else:  # var is missing
            # check if required inputs are available
            if item == "tdps" and "tdps_derived" not in df.columns:
                if _input_var_check(df, var1="tas", var2="hurs") == True:
                    print(
                        f"Calculating {item}_derived..."
                    )  # convert to logger when set-up
                    df["tdps_derived"] = _calc_dewpointtemp(df["tas"], df["hurs"])
                    # synergistic flag check
                    df = derive_synergistic_flag(df, "tdps_derived", "tas", "hurs")
            else:
                print("tdps_derived is present in station, no derivation necessary.")
                continue

            if item == "hurs" and _input_var_check(df, var1="tas", var2="tdps") == True:
                print(f"Calculating {item}_derived...")  # convert to logger when set-up
                df["hurs_derived"] = _calc_relhumid(df["tas"], df["tdps"])
                # synergistic flag check
                df = derive_synergistic_flag(df, "hurs_derived", "tas", "tdps")

            elif (
                item == "hurs"
                and _input_var_check(df, var1="tas", var2="tdps_derived") == True
            ):
                print(
                    f"Calculating {item}_derived ..."
                )  # convert to logger when set-up
                df["hurs_derived"] = _calc_relhumid(df["tas"], df["tdps_derived"])
                # synergistic flag check
                df = derive_synergistic_flag(df, "hurs_derived", "tas", "tdps_derived")

            elif (
                item == "tas" and _input_var_check(df, var1="hurs", var2="tdps") == True
            ):
                print(f"Calculating {item}_derived...")  # convert to logger when set-up
                df["tas_derived"] = _calc_airtemp(df["hurs"], df["tdps"])
                # synergistic flag check
                df = derive_synergistic_flag(df, "tas_derived", "hurs", "tdps")

            elif (
                item == "tas"
                and _input_var_check(df, var1="hurs", var2="tdps_derived") == True
            ):
                print(
                    f"Calculating {item}_derived ...."
                )  # convert to logger when set-up
                df["tas_derived"] = _calc_airtemp(df["hurs"], df["tdps_derived"])
                # synergistic flag check
                df = derive_synergistic_flag(df, "tas_derived", "tas", "tdps_derived")
            else:
                print(
                    f"{item} is missing the required input variables. {item}_derived not calculated."
                )  # convert to logger when set-up

        # TODO: attribute modification to denote it was derived

    return df


def _input_var_check(df: pd.DataFrame, var1: str, var2: str) -> bool:
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

    Returns
    -------
    bool
        True if all required input vars present; False if not
    """

    if var1 in df.columns and var2 in df.columns:
        return True
    else:
        return False


def derive_synergistic_flag(
    df: pd.DataFrame, var_to_flag: str, var1: str, var2: str
) -> pd.DataFrame:
    """Synergistically flags the derived variable if the input variables also have flags.

    Parameters
    ----------
    df : pd.DataFrame
        input df to identify flags
    var_to_flag : str
        name of variable to check and flag
    var1 : str
        name of secondary input var 1
    var2 : str
        name of secondary input var 2

    Returns
    -------
    df : pd.DataFrame
        df with synergistic flags applied, if applicable
    """
    # set up _eraqc variable for new derived variable
    df[var_to_flag + "_eraqc"] = np.nan

    # identify if var 1 has flags
    if len(df[var1 + "_eraqc"].unique()) > 1:
        # flags are present
        df.loc[df[var1 + "_eraqc"] > 0, var_to_flag + "_eraqc"] = (
            38  # see qaqc flag meanings
        )

    if len(df[var2 + "_eraqc"].unique()) > 1:
        df.loc[df[var2 + "_eraqc"] > 0, var_to_flag + "_eraqc"] = (
            38  # see qaqc flag meanings
        )

    return df


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
    return np.round(tdps, decimals=3)


def _calc_airtemp(hurs: pd.Series, tdps: pd.Series) -> pd.Series:
    """Calculate air temperature

    Parameters
    ----------
    hurs : pd.Series
        relative humidity, % or 0-100
    tdps : pd.Series
        dewpoint temperature, K

    Returns
    -------
    tas : pd.Series
        air temperature, K

    Notes
    ------
    [1] August-Roche-Magnus Approximation
    """

    # tdps must be in degC, not K for this equation
    tdps_degC = tdps - 273.15

    # apply approximation to calculate tas in degC
    tas_degC = (
        243.04
        * (((17.625 * tdps_degC) / (243.04 + tdps_degC)) - np.log(hurs / 100))
        / (17.625 + np.log(hurs / 100) - ((17.625 * tdps_degC) / (243.04 + tdps_degC)))
    )

    # convert back to K
    tas_K = tas_degC + 273.15

    return np.round(tas_K, decimals=3)


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
    return np.round(hurs, decimals=3)
