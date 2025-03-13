"""
This is a script where Stage 3: QA/QC function(s) whole station checks. 
For use within the PIR-19-006 Historical Obsevations Platform.
"""

## Import Libraries
import boto3
import geopandas as gp
import shapely
import numpy as np
import pandas as pd
import shapely
import scipy.stats as stats
import urllib
import requests

# New logger function
from log_config import logger

try:
    from qaqc_utils import *
except Exception as e:
    logger.debug("Error importing qaqc_utils: {}".format(e))

try:
    from qaqc_plot import flagged_timeseries_plot
except Exception as e:
    logger.debug("Error importing flagged_timeseries_plot: {}".format(e))

# if __name__ == "__main__":
wecc_terr = (
    "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
)
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"


# ======================================================================
## Part 1a functions (whole station/network)
## Note: QA/QC functions in part 1a of whole station checks do not proceed through QA/QC if failure occurs


# ----------------------------------------------------------------------
# missing value check: double check that all missing value observations are converted to NA before QA/QC
def qaqc_missing_vals(df, verbose=False):
    """
    Test for any errant missing values that made it through cleaning and converts missing values to NaNs.
    Searches for missing values in 3_qaqc_data/missing_data_flags.csv.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    verbose : bool, optional
        if True, returns runtime output to terminal

    Returns
    -------
    If QAQC is successful, returns a dataframe with flagged values
    If QAQC fails, returns None
    """

    logger.info("Running: qaqc_missing_vals")

    missing_vals = pd.read_csv("missing_data_flags.csv")

    vars_to_remove = ["qc", "duration", "method"]
    all_vars = [
        var
        for var in df.columns
        if var
        not in [
            "lon",
            "lat",
            "time",
            "day",
            "hour",
            "month",
            "year",
            "date",
            "elevation",
            "station",
            "anemometer_height_m",
            "thermometer_height_m",
        ]
    ]
    obs_vars = [
        var
        for var in all_vars
        if not any(True for item in vars_to_remove if item in var)
    ]

    try:
        for item in obs_vars:
            # pull missing values which are appropriate for the range of real values for each variable
            missing_codes = missing_vals.loc[
                missing_vals["variable"].str.contains(item)
                | missing_vals["variable"].str.contains("all")
            ]

            # values in column that == missing_flag values, replace with NAs
            # note numerical vals converted to strings first to match missing_flag formatting
            df[item] = np.where(
                df[item].astype(str).isin(missing_codes["missing_flag"]),
                float("NaN"),
                df[item],
            )

            logger.info(
                "Updating missing values for: {}".format(item),
            )
            return df

    except Exception as e:
        logger.info(e)
        return None


# ----------------------------------------------------------------------
# missing spatial coords (lat-lon)
def qaqc_missing_latlon(df, verbose=False):
    """
    Test for missing latitude / longitude values for a station.
    Checks if latitude and longitude is missing for a station.
    If missing, station is flagged to not proceed through QA/QC.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    verbose : bool, optional
        if True, returns runtime output to terminal

    Returns
    -------
    If QAQC is successful, returns a dataframe with flagged values
    If QAQC fails, returns None
    """
    logger.info("Running: qaqc_missing_latlon")

    # latitude or longitude
    variables = list(df.columns)
    if "lon" not in variables or "lat" not in variables:
        return None

    if df["lat"].isnull().all():
        return None

    if df["lon"].isnull().all():
        return None

    return df


# ----------------------------------------------------------------------
# in bounds of WECC
def qaqc_within_wecc(df, verbose=False):
    """
    Test for whether station is within terrestrial & marine WECC boundaries.
    If outside of boundaries, station is flagged to not proceed through QA/QC.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    verbose : bool, optional
        if True, returns runtime output to terminal

    Returns
    -------
    If QAQC is successful, returns a dataframe with flagged values
    If QAQC fails, returns None
    """
    logger.info("Running: qaqc_within_wecc")

    t = gp.read_file(wecc_terr).iloc[0].geometry  ## Read in terrestrial WECC shapefile.
    m = gp.read_file(wecc_mar).iloc[0].geometry  ## Read in marine WECC shapefile.
    pxy = shapely.geometry.Point(df["lon"].mean(), df["lat"].mean())
    if pxy.within(t) or pxy.within(m):
        return df
    else:
        return None


# ----------------------------------------------------------------------
# elevation
def _grab_dem_elev_m(lats_to_check, lons_to_check, verbose=False):
    """
    If elevation is missing, retrieves elevation value from the USGS Elevation Point Query Service,
    lat lon must be in decimal degrees (which it is after cleaning).

    Parameters
    ----------
    lats_to_check : list of float
        list of latitude values to input for DEM elevation retrieval
    lons_to_check : list of float
        list of longitude values to input for DEM elevation retrieval
    verbose : bool, optional
        if True, returns runtime output to terminal

    Returns
    -------
    dem_elev_short : float
        elevation infill value

    References
    ----------
    [1] Modified from: https://gis.stackexchange.com/questions/338392/getting-elevation-for-multiple-lat-long-coordinates-in-python
    """

    url = r"https://epqs.nationalmap.gov/v1/json?"

    dem_elev_short = np.ones_like(lats_to_check) * np.nan

    for i, lat, lon in zip(range(len(lats_to_check)), lats_to_check, lons_to_check):
        # define rest query params
        params = {"output": "json", "x": lon, "y": lat, "units": "Meters"}
        # format query string and return value
        result = requests.get((url + urllib.parse.urlencode(params)))
        dem_elev_long = float(result.json()["value"])
        # make sure to round off lat-lon values so they are not improbably precise for our needs
        dem_elev_short[i] = np.round(dem_elev_long, decimals=2)

    return dem_elev_short.astype("float")


# ----------------------------------------------------------------------
def qaqc_elev_infill(df, verbose=False):
    """Test if elevation is NA/missing.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    verbose : bool, optional
        if True, returns runtime output to terminal

    Returns
    -------
    If QAQC is successful, returns a dataframe with flagged values (see below for flag meaning)
    If QAQC fails, returns None

    Notes
    -----
    1. Three in-filling scenarios
        - If all values are nan, in-fill from DEM
        - If some values are missing, in-fill from station for consistency (e.g., OtherISD, ASOSAWOS)
        - If all values are nan, and DEM fails to retrieve (over ocean), set to 0 (e.g., NDBC, MARITIME)
    2. Elevation qa/qc flags: observations with these flags should continue through QA/QC
    Flag meaning : 3,qaqc_elev_infill,Elevation infilled from DEM (USGS 3DEP)
    Flag meaning : 4,qaqc_elev_infill,Elevation infilled from station
    Flag meaning : 5,qaqc_elev_infill,Elevation manually infilled for buoy out-of-range of DEM
    """

    logger.info("Running: qaqc_elev_infill")

    # first check to see if any elev value is missing
    if df["elevation"].isnull().any() == True:
        # all elevation values are reported as nan (some ndbc/maritime)

        if df["elevation"].isnull().values.all() == True:
            try:  # in-fill if value is missing
                # check if lat-lon has changed over time
                nan_lats = df["lat"].unique()
                nan_lons = df["lon"].unique()

                if (len(nan_lats) == 1) and (
                    len(nan_lons) == 1
                ):  # single lat-lon pair for missing elevs
                    try:
                        dem_elev_value = _grab_dem_elev_m(
                            list(nan_lats), list(nan_lons)
                        )
                        df.loc[df["elevation"].isnull() == True, "elevation_eraqc"] = (
                            3  # see era_qaqc_flag_meanings.csv
                        )
                        df.loc[df["elevation"].isnull() == True, "elevation"] = float(
                            dem_elev_value
                        )

                    except:  # some buoys out of range of dem (past coastal range) report nan elevation, manually set to 0.00m and flag
                        df.loc[df["elevation"].isnull() == True, "elevation_eraqc"] = (
                            5  # see era_qaqc_flag_meanings.csv
                        )
                        df.loc[df["elevation"].isnull() == True, "elevation"] = float(
                            0.00
                        )  # manual infilling for buoys

                else:  # multiple pairs of lat-lon for missing elevs
                    for ilat in nan_lats:
                        for ilon in nan_lons:
                            dem_elev_value = _grab_dem_elev_m(ilat, ilon)
                            df.loc[
                                (df["lat"] == ilat) & (df["lon"] == ilon),
                                "elevation_eraqc",
                            ] = 3  # see era_qaqc_flag_meanings.csv
                            df.loc[
                                (df["lat"] == ilat) & (df["lon"] == ilon), "elevation"
                            ] = float(dem_elev_value)

            except:  # elevation cannot be obtained from DEM
                logger.info("Elevation cannot be in-filled")
                return None

        # some stations have a single/few nan reported, in-fill from station for consistency
        else:  # multiple values for elevation, infill each instance if missing/incorrectly coded (e.g., zeros when shouldnt be)
            try:
                # locate all instances of nan values as elevation codes
                nan_coded = df[df["elevation"].isnull() == True]
                nan_lats = nan_coded["lat"].unique()
                nan_lons = nan_coded["lon"].unique()

                if (len(nan_lats) == 1) and (
                    len(nan_lons) == 1
                ):  # single lat-lon pair for missing elevs
                    if (nan_lats[0] == df["lat"].iloc[0]) & (
                        nan_lons[0] == df["lon"].iloc[0]
                    ):  # single set of lat-lons matches station, infill from station
                        df.loc[df["elevation"].isnull() == True, "elevation_eraqc"] = (
                            4  # see era_qaqc_flag_meanings.csv
                        )
                        df.loc[df["elevation"].isnull() == True, "elevation"] = df[
                            "elevation"
                        ].iloc[0]
                    else:  # lat-lon of missing elev does not match station lat-lon (has shifted), infill from dem
                        dem_elev_value = _grab_dem_elev_m(nan_lats[0], nan_lons[0])
                        df.loc[df["elevation"].isnull() == True, "elevation_eraqc"] = (
                            3  # see era_qaqc_flag_meanings.csv
                        )
                        df.loc[df["elevation"].isnull() == True, "elevation"] = float(
                            dem_elev_value
                        )

                else:  # multiple pairs of lat-lon for missing elevs
                    for ilat in nan_lats:
                        for ilon in nan_lons:
                            dem_elev_value = _grab_dem_elev_m(ilat, ilon)
                            df.loc[
                                (df["lat"] == ilat) & (df["lon"] == ilon),
                                "elevation_eraqc",
                            ] = 3  # see era_qaqc_flag_meanings.csv
                            df.loc[
                                (df["lat"] == ilat) & (df["lon"] == ilon), "elevation"
                            ] = float(dem_elev_value)

            # elevation cannot be in-filled
            except:
                logger.info("Elevation cannot be in-filled")
                return None
    return df


# ----------------------------------------------------------------------
def qaqc_elev_range(df, verbose=False):
    """
    Checks if valid elevation value is outside of range of reasonable values for WECC region.
    If outside range, station is flagged to not proceed through QA/QC.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    verbose : bool, optional
        if True, returns runtime output to terminal

    Returns
    -------
    If QAQC is successful, returns a dataframe with flagged values
    If QAQC fails, returns None

    Notes
    -----
    1. Death Valley is 282 (85.9 m) below sea level
    2. Denali is ~6190 m above sea level
    3. Opting to build in a 10 m buffer for range check
    """

    logger.info("Running: qaqc_elev_range")

    # If value is present but outside of reasonable value range
    if (df["elevation"].values.any() < -95.0) or (
        df["elevation"].values.any() > 6210.0
    ):
        logger.info(
            "Station out of range for WECC -- station does not proceed through QAQC"
        )
        return None

    # Elevation value is present and within reasonable value range
    else:
        df = df
        logger.info(
            "Elevation values post-infilling/correcting: {}".format(
                df["elevation"].unique()
            )
        )

    return df


# ======================================================================
## Part 1b functions (whole station/network)
## Note: QA/QC functions in part 1b of whole station checks proceed through QA/QC if failure occurs


# ----------------------------------------------------------------------
## sensor height - air temperature
def qaqc_sensor_height_t(df, verbose=False):
    """
    Checks if temperature sensor height is within 2 meters above surface +/- 1/3 meter tolerance.
    If missing or outside range, temperature value for station is flagged to not proceed through QA/QC.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    verbose : bool, optional
        if True, returns runtime output to terminal

    Returns
    -------
    If QAQC is successful, returns a dataframe with flagged values (see below for flag meaning)
    If QAQC fails, returns None

    Notes
    ------
    1. qaqc_sensor_height_t function moved into v2 of this data product, as many networks do not report sensor height, leaving many stations excluded from this check
    Flag meaning : 6,qaqc_sensor_height_t,Thermometer height missing
    Flag meaning : 7,qaqc_sensor_height_t,Thermometer height not 2 meters
    """

    logger.info("Running: qaqc_sensor_height_t")

    try:
        # Check if thermometer height is missing
        isHeightMissing = df["thermometer_height_m"].isnull().any()

        if isHeightMissing:
            df["tas_eraqc"] = 6  # see era_qaqc_flag_meanings.csv
            logger.info(
                "Thermometer height is missing -- air temperature will be excluded from all QA/QC checks",
            )
        else:
            isHeightWithin = np.logical_and(
                df["thermometer_height_m"] >= (2 - 1 / 3),
                df["thermometer_height_m"] <= (2 + 1 / 3),
            )

            # Thermometer height present but outside 10m +/- tolerance
            if not isHeightWithin.all():
                # df.loc[:, 'tas_eraqc'] = 7 # see era_qaqc_flag_meanings.csv
                df["tas_eraqc"] = 7  # see era_qaqc_flag_meanings.csv
                logger.info(
                    "Thermometer height is not 2 m -- air temperature will be excluded from all QA/QC checks",
                )

        return df
    except Exception as e:
        logger.info(
            "qaqc_sensor_height_t failed with Exception: {}".format(e),
        )
        return None


# ----------------------------------------------------------------------
## sensor height - wind
## NOTE: qaqc_sensor_height_w function moved into v2 of this data product, as many networks do not report sensor height, leaving many stations excluded from this check
def qaqc_sensor_height_w(df, verbose=False):
    """
    Checks if wind sensor height is within 10 meters above surface +/- 1/3 meter tolerance.
    If missing or outside range, wind speed and direction values for station are flagged to not proceed through QA/QC.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    verbose : bool, optional
        if True, returns runtime output to terminal

    Returns
    -------
    If QAQC is successful, returns a dataframe with flagged values (see below for flag meaning)
    If QAQC fails, returns None

    Notes
    ------
    1. qaqc_sensor_height_w function moved into v2 of this data product, as many networks do not report sensor height, leaving many stations excluded from this check
    Flag meaning : 8,qaqc_sensor_height_w,Anemometer height missing
    Flag meaning : 9,qaqc_sensor_height_w,Anemometer height not 10 meters
    """

    logger.info("Running: qaqc_sensor_height_w")

    try:
        # Check if anemometer height is missing
        isHeightMissing = df["anemometer_height_m"].isnull().any()

        if isHeightMissing:
            # df.loc[:,'sfcWind_eraqc'] = 8 # see era_qaqc_flag_meanings.csv
            # df.loc[:,'sfcWind_dir_eraqc'] = 8
            df["sfcWind_eraqc"] = 8  # see era_qaqc_flag_meanings.csv
            df["sfcWind_dir_eraqc"] = 8
            logger.info(
                "Anemometer height is missing -- wind speed and direction will be excluded from all QA/QC checks",
            )

        else:  # sensor height present
            # Check if anemometer height is within 10 m +/- 1/3 m
            isHeightWithin = df["anemometer_height_m"][0] >= (10 - 1 / 3) and df[
                "anemometer_height_m"
            ][0] <= (10 + 1 / 3)
            # Anemometer height present but outside 10m +/- tolerance
            if not isHeightWithin.all():
                df["sfcWind_eraqc"] = 9  # see era_qaqc_flag_meanings.csv
                df["sfcWind_dir_eraqc"] = 9
                logger.info(
                    "Anemometer height is not 10 m -- wind speed and direction will be excluded from all QA/QC checks",
                )
        return df

    except Exception as e:
        logger.info(
            "qaqc_sensor_height_w failed with Exception: {}".format(e),
        )
        return None


# ----------------------------------------------------------------------
## flag values outside world records for North America
def qaqc_world_record(df, verbose=False):
    """
    Checks if variables are outside North American world records.
    If outside minimum or maximum records, flags values.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    verbose : bool, optional
        if True, returns runtime output to terminal

    Returns
    -------
    If QAQC is successful, returns a dataframe with flagged values (see below for flag meaning)
    If QAQC fails, returns None

    Notes
    ------
    Flag meaning : 11,qaqc_world_record,Value outside of world record range

    References
    ----------
    [1] World records from HadISD protocol, cross-checked with WMO database
    [2] https://wmo.asu.edu/content/world-meteorological-organization-global-weather-climate-extremes-archive
    [3] Solar radiation specific: Rupp et al. 2022, Slater 2016
    [4] https://www.ncei.noaa.gov/access/monitoring/scec/records
    [5] https://www.weather.gov/media/owp/oh/hdsc/docs/TP2.pdf
    """

    logger.info("Running: qaqc_world_record")

    try:
        T_X = {"North_America": 329.92}  # temperature, K
        T_N = {"North_America": 210.15}  # temperature, K
        D_X = {"North_America": 329.85}  # dewpoint temperature, K
        D_N = {"North_America": 173.15}  # dewpoint temperature, K
        W_X = {"North_America": 113.2}  # wind speed, m/s
        W_N = {"North_America": 0.0}  # wind speed, m/s
        R_X = {"North_America": 1500}  # solar radiation, W/m2
        R_N = {"North_America": -5}  # solar radiation, W/m2

        # for other non-record variables (wind direction, humidity)
        N_X = {"North_America": 360}  # wind direction, degrees
        N_N = {"North_America": 0}  # wind direction, degrees
        H_X = {"North_America": 100}  # humidity, max
        H_N = {"North_America": 0}  # humidity, min
        E_X = {"North_America": 6210.0}  # elevation, m
        E_N = {"North_America": -100}  # elevation, m

        # pressure, with elevation options
        S_X = {"North_America": 108330}  # pressure, Pa
        S_N = {"North_America": 87000}  # sea level pressure only, Pa
        SALT_N = {
            "North_America": 45960
        }  # non-sea level pressure, Pa, reduced min based on max elevation (6190 m)

        # precipitation, with variations depending on reporting interval
        P_X = {"North_America": 656}  # precipitation, mm, 24-hr rainfall
        PALT5_X = {"North_America": 31.8}  # precipitation, mm, 5-min rainfall
        PALT15_X = {"North_America": 25.4}  # precipitation, mm, 15-min rainfall
        P_N = {"North_America": 0}  # precipitaiton, mm

        maxes = {
            "tas": T_X,
            "tdps": D_X,
            "tdps_derived": D_X,
            "sfcWind": W_X,
            "sfcWind_dir": N_X,
            "psl": S_X,
            "ps": S_X,
            "ps_derived": S_X,
            "ps_altimeter": S_X,
            "rsds": R_X,
            "pr": P_X,
            "pr_5min": PALT5_X,
            "pr_15min": PALT15_X,
            "pr_1h": P_X,
            "pr_24h": P_X,
            "pr_localmid": P_X,
            "hurs": H_X,
            "elevation": E_X,
        }
        mins = {
            "tas": T_N,
            "tdps": D_N,
            "tdps_derived": D_N,
            "sfcWind": W_N,
            "sfcWind_dir": N_N,
            "psl": S_N,
            "ps": SALT_N,
            "ps_derived": SALT_N,
            "ps_altimeter": SALT_N,
            "rsds": R_N,
            "pr": P_N,
            "pr_5min": P_N,
            "pr_15min": P_N,
            "pr_1h": P_N,
            "pr_24h": P_N,
            "pr_localmid": P_N,
            "hurs": H_N,
            "elevation": E_N,
        }

        # variable names to check against world record limits
        wr_vars = [
            "tas",
            "tdps",
            "tdps_derived",
            "sfcWind",
            "sfcWind_dir",
            "ps",
            "psl",
            "ps_altimeter",
            "ps_derived",
            "rsds",
            "pr",
            "pr_5min",
            "pr_15min",
            "pr_1h",
            "pr_24h",
            "pr_24h",
            "hurs",
            "elevation",
        ]
        for var in wr_vars:
            if var in list(df.columns):
                df_valid = grab_valid_obs(df, var)  # subset for valid obs
                isOffRecord = np.logical_or(
                    df_valid[var] < mins[var]["North_America"],
                    df_valid[var] > maxes[var]["North_America"],
                )
                if isOffRecord.any():
                    df.loc[isOffRecord, var + "_eraqc"] = (
                        11  # see era_qaqc_flag_meanings.csv
                    )
                    logger.info(
                        "Flagging {} observations exceeding world/regional records: {}".format(
                            len(isOffRecord), var
                        )
                    )

        return df
    except Exception as e:
        logger.info(
            "qaqc_world_record failed with Exception: {}".format(e),
        )
        return None


# ----------------------------------------------------------------------
## final summary stats of flagged variables and percentage of coverage
def flag_summary(df, verbose=False, local=False):
    """Generates summary of flags set on all QAQC tests.
    Returns list of unique flag values for each variable
    Returns % of total obs per variable that was flagged
    Information is included in log_file.

    Also produces figures of full timeseries

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    verbose : bool, optional
        if True, returns runtime output to terminal
    local : bool, optional
        if True, saves plots to local folder in addition to AWS

    Returns
    -------
    None
        This function does not return a value
    """

    logger.info("Running: flag_summary")

    # identify _eraqc variables
    eraqc_vars = [var for var in df.columns if "_eraqc" in var]
    obs_vars = [item.split("_e")[0] for item in eraqc_vars]

    for var in eraqc_vars:
        logger.info(
            "Flags set on {}: {}".format(var, df[var].unique()),
        )  # unique flag values
        logger.info(
            "Coverage of {} obs flagged: {}% of obs".format(
                var,
                round((len(df.loc[(df[var].isnull() == False)]) / len(df)) * 100, 2),
            )
        )  # % of coverage flagged

    for var in obs_vars:
        try:
            flagged_timeseries_plot(df, var)
        except Exception as e:
            logger.info(
                "flagged_timeseries_plot failed for {} with Exception: {}".format(
                    var, e
                )
            )
