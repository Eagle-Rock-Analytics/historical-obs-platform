"""merge_derive_missing.py

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
import logging
import inspect


## Identify vars that can be derived
def merge_derive_missing_vars(
    df: pd.DataFrame, var_attrs: dict, logger: logging.Logger
) -> tuple[pd.DataFrame, dict] | None:
    """
    Identifies if any variables can be derived with other input variables.
    If success, variable is derived in the correct unit, attribtues are updated,
    and any flags from the input variables are synergistically flagged.
    If failure, variable is not derived.

    Parameters
    ----------
    df : pd.DataFrame
        input dataframe
    var_attrs : dict
        attributes of input variables

    Returns
    -------
    If success: pd.DataFrame with newly added derived variable, and updated variable attributes
    If failure: None
    """
    logger.info(f"{inspect.currentframe().f_code.co_name}: Starting...")

    # vars that can be derived
    derive_vars = [
        "hurs",
        "tas",
    ]  # tdps not included here, kept separate because of tdps_derived

    # initialize update vars dictionary
    new_var_attrs = var_attrs.copy()

    try:
        # var is missing
        # check if required inputs are available
        if "tdps" not in df.columns and "tdps_derived" not in df.columns:
            if _input_var_check(df, var1="tas", var2="hurs") == True:
                logger.info("Calculating tdps_derived...")  # convert to logger when set-up
                df["tdps_derived"] = _calc_dewpointtemp(df["tas"], df["hurs"])
                # synergistic flag check
                df = derive_synergistic_flag(df, "tdps_derived", "tas", "hurs")
                # add new variable attributes
                new_var_attrs = _add_derived_var_attrs(
                    derived_var="tdps_derived",
                    source_var="tdps",
                    input_vars=["tas", "hurs"],
                    var_attrs=new_var_attrs,
                )
                logger.info("Successfully calculated tdps_derived")

        else:
            logger.info("tdps_derived is present in station, no derivation necessary.")

        # first check if station has any vars that can be derived
        for item in derive_vars:
            if item in df.columns:
                logger.info(
                    f"{item} is present in station, no derivation necessary."
                )  # convert to logger when set-up
                continue

            if item == "hurs" and _input_var_check(df, var1="tas", var2="tdps") == True:
                logger.info(f"Calculating {item}_derived...")  # convert to logger when set-up
                df["hurs_derived"] = _calc_relhumid(df["tas"], df["tdps"])
                # synergistic flag check
                df = derive_synergistic_flag(df, "hurs_derived", "tas", "tdps")
                # add new variable attributes
                new_var_attrs = _add_derived_var_attrs(
                    derived_var="hurs_derived",
                    source_var="hurs",
                    input_vars=["tas", "tdps"],
                    var_attrs=new_var_attrs,
                )
                logger.info(f"Successfully calculated {item}_derived")

            elif (
                item == "hurs"
                and _input_var_check(df, var1="tas", var2="tdps_derived") == True
            ):
                logger.info(
                    f"Calculating {item}_derived ..."
                )  # convert to logger when set-up
                df["hurs_derived"] = _calc_relhumid(df["tas"], df["tdps_derived"])
                # synergistic flag check
                df = derive_synergistic_flag(df, "hurs_derived", "tas", "tdps_derived")
                # add new variable attributes
                new_var_attrs = _add_derived_var_attrs(
                    derived_var="tdps_derived",
                    source_var="tdps",
                    input_vars=["tas", "hurs"],
                    var_attrs=new_var_attrs,
                )
                logger.info(f"Successfully calculated {item}_derived")

            elif (
                item == "tas" and _input_var_check(df, var1="hurs", var2="tdps") == True
            ):
                logger.info(f"Calculating {item}_derived...")  # convert to logger when set-up
                df["tas_derived"] = _calc_airtemp(df["hurs"], df["tdps"])
                # synergistic flag check
                df = derive_synergistic_flag(df, "tas_derived", "hurs", "tdps")
                # add new variable attributes
                new_var_attrs = _add_derived_var_attrs(
                    derived_var="tas_derived",
                    source_var="tas",
                    input_vars=["hurs", "tdps"],
                    var_attrs=new_var_attrs,
                )
                logger.info(f"Successfully calculated {item}_derived")

            elif (
                item == "tas"
                and _input_var_check(df, var1="hurs", var2="tdps_derived") == True
            ):
                logger.info(
                    f"Calculating {item}_derived ...."
                )  # convert to logger when set-up
                df["tas_derived"] = _calc_airtemp(df["hurs"], df["tdps_derived"])
                # synergistic flag check
                df = derive_synergistic_flag(df, "tas_derived", "tas", "tdps_derived")
                # add new variable attributes
                new_var_attrs = _add_derived_var_attrs(
                    derived_var="tas_derived",
                    source_var="tas",
                    input_vars=["hurs", "tdps_derived"],
                    var_attrs=new_var_attrs,
                )
                logger.info(f"Successfully calculated {item}_derived")

            else:
                logger.info(
                    f"{item} is missing the required input variables. {item}_derived not calculated."
                )  # convert to logger when set-up

        logger.info(f"{inspect.currentframe().f_code.co_name}: Completed successfully")
        return df, new_var_attrs

    except Exception as e:
        logger.error(f"{inspect.currentframe().f_code.co_name}: Failed")
        raise e

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

    Notes
    -----
    Flag meaning : 38,derive_synergistic_flag,At least one input variable to derived variable has a flag placed. Input variable and derived variable are synergistically flagged

    """
    # set up _eraqc variable for new derived variable
    df[var_to_flag + "_eraqc"] = np.nan

    # Convert string qaqc column to float
    # This column shouldn't have been strings in the first place, but oh well.
    for var in [var1, var2]:
        col = var + "_eraqc"
        df[col] = df[col].replace(
            ["nan", ""], np.nan
        )  # Replace string "nan" and empty string "" to np.nan
        df[col] = pd.to_numeric(
            df[col], errors="coerce"
        )  # Convert string integers ("28") to actual integers (28)

    # identify if var 1 has flags
    if len(df[var1 + "_eraqc"].unique()) > 1:
        # flags are present
        df.loc[df[var1 + "_eraqc"] > 0, var_to_flag + "_eraqc"] = (
            38  # see 3_qaqc_data/era_qaqc_flag_meanings.csv
        )

    if len(df[var2 + "_eraqc"].unique()) > 1:
        df.loc[df[var2 + "_eraqc"] > 0, var_to_flag + "_eraqc"] = (
            38  # see 3_qaqc_data/era_qaqc_flag_meanings.csv
        )

    return df


def _add_derived_var_attrs(
    derived_var: str, source_var: str, input_vars: list[str], var_attrs: dict
) -> dict:
    """Creates data attributes for new derived variable and adds to var_attrs.

    Parameters
    ----------
    derived_var : str
        variable name of new derived variable
    source_var : str
        variable name of the variable it "derives"
    input_vars : list[str]
        variable names of input variable
    var_attrs : dict
        attributes for all variables

    Returns
    -------
    var_attrs : dict
        updated variable attributes dictionary with new vars
    """

    # support for naming, units
    if source_var == "tdps":
        long_name = "derived_dew_point_temperature"
        units = "K"
    elif source_var == "tas":
        long_name = "derived_air_temperature"
        units = "K"
    elif source_var == "hurs":
        long_name = "derived_relative_humidity"
        units = "percent"

    # add new attributes -- var_attrs are stored as dict of each var dict
    derived_var_dict = {
        "long_name": long_name,
        "units": units,
        "ancillary_variables": f"{input_vars[0]}, {input_vars[1]}",
        "comment": "Derived in merge_derive_missing_vars.",
    }

    # add new var dictionary to existing var_attrs dict
    var_attrs[derived_var] = derived_var_dict
    return var_attrs


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

    Notes
    -----
    Rounded to 3 decimal places to be consistent with input raw data sig figs
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

    References
    ----------
    [1] August-Roche-Magnus Approximation

    Notes
    -----
    Rounded to 3 decimal places to be consistent with input raw data sig figs
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

    Notes
    -----
    Rounded to 3 decimal places to be consistent with input raw data sig figs
    """

    es = 0.611 * np.exp(
        5423 * ((1 / 273) - (1 / tas))
    )  # calculates saturation vapor pressure using air temp
    e_vap = 0.611 * np.exp(
        5423 * ((1 / 273) - (1 / tdps))
    )  # calculates vapor pressure using dew point temp
    hurs = 100 * (e_vap / es)
    return np.round(hurs, decimals=3)
