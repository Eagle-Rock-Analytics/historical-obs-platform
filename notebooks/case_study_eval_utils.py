"""
case_study_eval_utils.py

Functions used to support case study analysis of extreme events.

Functions
---------
- known_issue_check: Identifies if station under evaluation has a known network issue.
- subset_eval_stns: Identifies stations to evaluate for specific V1 case study events.
- id_all_flags: Prints all unique values of all eraqaqc flags.
- event_info: Utility function to return useful information for a specific designated event.
- event_subset: Subsets for the event itself + buffer around to identify event.
- flags_during_event: Provides info on which flags were placed during event for evaluation.
- find_other_events: Event finder not tied to specified case study events.
- latlon_to_mercator_cartopy: Converts lat/lon coordinates to mercator for plotting.
- stn_visualize: Produces simple map of station relevant to event boundary.
- event_plot: Produces timeseries of variables that have flags placed.

Intended Use
------------
Used in the case study analysis for ease of comparison and tidy notebooks.
"""

from pyproj import CRS, Transformer
import geopandas as gpd
from geopandas import GeoDataFrame
from shapely.geometry import Point, Polygon
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.feature as cf
from matplotlib.ticker import MaxNLocator
import cartopy.crs as ccrs
import datetime
import sys
import os

# Import qaqc stage plot functions
sys.path.append(os.path.abspath("../scripts/3_qaqc_data"))
from qaqc_plot import flagged_timeseries_plot, _plot_format_helper, id_flag

CENSUS_SHP = "s3://wecc-historical-wx/0_maps/ca_counties/CA_Counties.shp"


def known_issue_check(network: str, var: str, stn: str):
    """
    Identifies if station under evaluation has a known network issue.
    At present, only prints out a statement if there is an issue.
    Eventually may want to do <something>

    Parameters
    ----------
    network : str
        Name of network to check
    var : str
        Name of variable to check
    stn : str
        Name of station to check

    Returns
    -------
    None

    Notes
    -----
    1. See "Known Network Issues for QA/QC Validation" planning doc.
    """
    print("Checking for known station issues...")

    # RAWS
    if network == "RAWS":
        if var == "tas":
            print(
                f"Known network issue for {network} {var}: values may be too high (on order of 10°F) if sun is shining strongly and winds are light."
            )

        elif var == "pr":
            print(
                f"Known network issue for {network} {var}: stations are not maintained in winter, instrument may freeze. Consider subsetting for May-September."
            )
            # V2 note: exclude RAWS data during specific notes -- would require new function to flag

    # SNOTEL
    if network == "SNOTEL":
        if var == "tas":
            print(
                f"Known network issue for {network} {var}: values may remain at exactly 0.0°C for two or more consecutive days. Should be caught by unusual_streaks."
            )
            print(
                f"Known network issue for {network} {var}: SNOTEL temperature sensors transition between mid-1990s and mid-2000s to new sensory type produces warm bias at \
            colder temperatures. Min temperature may be too high, max temperature may be too low."
            )
            # V2 note: trend analysis may identify these issues, nearest neighbor check could identify

    # ASOSAWOS + OtherISD
    if network == "ASOSAWOS":
        if var == "tdps":
            print(
                f"Known network issue for {network} {var}: values may be stuck at around 0.0°C, or have excessive mirror contamination. Should be caught by unsusual_streaks."
            )

    if network == "ASOSAWOS" or network == "OtherISD":
        if var == "pr":
            print(
                f"Known network issue for {network} {var}: ASOS network began installation in 1996, with poor instrumentation for measuring snowfall. Precipitation between \
            1980-1996 may be more likely to be flagged."
            )

    # CIMIS
    if network == "CIMIS":
        if var == "pr":
            print(
                f"Known network issue for {network} {var}: stations located in flat agricultural areas, sensor may be detecting sprinkler irrigation events. \
            Network does have stringent QC protocol."
            )
            # V2 note: nearest neighbor check could confirm

    # NDBC / MARITIME
    if network == "NDBC" or network == "MARITIME":
        print(
            f"Known network issue for {network}: some buoys have data past their known disestablishment dates. Should be caught by spurious_buoy_check."
        )

        if stn == "NDBC_46044":
            print(
                "Known network issue for NDBC_46044: buoy went adrift during reporting period. Confirm if data was flagged by QA/QC."
            )
            # V2 note: if not flagged, needs to be -- would require new function

        if (
            stn == "MARITIME_MTYC1"
            or stn == "MARITIME_MEYC1"
            or stn == "MARITIME_SMOC1"
            or stn == "MARITIME_ICAC1"
        ):
            print(
                f"Known network issue for {network} station {stn}: buoy was renamed and/or relocated. May cause issue for station proximity tests."
            )
            # V2 note: noted in qaqc_buoy_check but not handled -- would require new function

    return None


def subset_eval_stns(
    event_to_eval: str,
    stn_list: pd.DataFrame,
    specific_station: str | None = None,
    subset: int | None = None,
    return_stn_ids: bool = False,
) -> pd.DataFrame:
    """
    Identifies stations to evaluate for specific V1 QA/QC events.
    Option to subset to a more manageable number of random stations for initial evaluation.

    Parameters
    ----------
    event_to_eval : str
        options: santa_ana_wind, winter_storm, AR, mudslide, aug2020_heatwave, sep2020_heatwave, aug2022_heatwave, offshore_wind
    stn_list : pd.DataFrame
        station list
    specific_station : str, optional
        name of specific station to check
    subset : int, optional
        optional value to specify number of stations to return, useful for big events
    return_stn_ids : bool, optional
        Option to return string names for ease of use

    Returns
    -------
    eval_stns : pd.DataFrame
        Subset of stations relevant to a desired extreme event for analysis

    To Dos
    -------
    1. Validation check on event_to_eval options
    2. Need an option for "WECC wide" (or no spatial subsetting)
    3. Update station list being used here -- SNOTEL dates have been fixed.
    """

    event_flags = []
    event_flags.append("all")
    event_flags.append(event_to_eval)

    # grab stations per event
    event_stns = stn_list[stn_list["event_type"].isin(event_flags)]

    # exclude "manual check on end date" for time being -- SNOTEl stations all have 2100 as their end date regardless of when data actually ends
    mask = event_stns["notes"] == "manual check on end date"
    event_stns = event_stns[~mask]
    # print('{} potential stations available for evaluation for {} event!'.format(len(event_stns), event_to_eval))

    # identify stations in geographic region we are looking for
    ca_county = gpd.read_file(CENSUS_SHP)

    # different areas based on events
    if event_to_eval == "santa_ana_wind":
        counties_to_grab = [
            "Los Angeles",
            "Orange",
            "San Diego",
            "San Bernardino",
            "Riverside",
            "Ventura",
            "Santa Barbara",
        ]

    elif event_to_eval == "winter_storm":
        # focus on Northern/Central/Bay Area to begin with // WECC wide
        counties_to_grab = [
            "Butte",
            "Colusa",
            "Del Norte",
            "Glenn",
            "Humboldt",
            "Lake",
            "Lassen",
            "Mendocino",
            "Modoc",
            "Nevada",
            "Plumas",
            "Shasta",
            "Sierra",
            "Siskiyou",
            "Tehama",
            "Trinity",
            "Alpine",
            "Amador",
            "Calaveras",
            "El Dorado",
            "Fresno",
            "Inyo",
            "Kings",
            "Madera",
            "Mariposa",
            "Merced",
            "Mono",
            "Placer",
            "Sacramento",
            "San Joaquin",
            "Stanislaus",
            "Sutter",
            "Yuba",
            "Tulare",
            "Tuolumne",
            "Yolo",
            "Alameda",
            "Contra Costa",
            "Marin",
            "Monterey",
            "Napa",
            "San Benito",
            "San Francisco",
            "San Mateo",
            "Santa Clara",
            "Santa Cruz",
            "Solano",
            "Sonoma",
        ]

    elif event_to_eval == "mudslide":
        counties_to_grab = ["Santa Barbara"]

    elif event_to_eval == "AR":
        counties_to_grab = []  # CA

    elif event_to_eval == "aug2020_heatwave":
        counties_to_grab = []  # CA

    elif event_to_eval == "sep2020_heatwave":
        counties_to_grab = [
            "San Luis Obispo",
            "Kern",
            "San Bernadino",
            "Santa Barbara",
            "Ventura",
            "Los Angeles",
            "Orange",
            "Riverside",
            "San Diego",
            "Imperial",
        ]

    elif event_to_eval == "aug2022_heatwave":
        # August 2022 -- Labor Day Heatwave "aug2022_heatwave"
        counties_to_grab = [
            "San Luis Obispo",
            "Kern",
            "San Bernadino",
            "Santa Barbara",
            "Ventura",
            "Los Angeles",
            "Orange",
            "Riverside",
            "San Diego",
            "Imperial",
        ]

    elif event_to_eval == "offshore_wind":
        counties_to_grab = [
            "San Diego",
            "Orange",
            "Los Angeles",
            "Ventura",
            "Santa Barbara",
            "San Luis Obispo",
            "Monterey",
            "Santa Cruz",
            "San Mateo",
            "Santa Clara",
            "Alameda",
            "San Francisco",
            "Contra Costa",
            "Solano",
            "Marin",
            "Sonoma",
            "Mendocino",
            "Humboldt",
            "Del Norte",
        ]

    target_counties = ca_county[ca_county["NAME"].isin(counties_to_grab)]
    target_counties = GeoDataFrame(target_counties, geometry=target_counties.geometry)

    geometry = [
        Point(latlon_to_mercator_cartopy(lat, lon))
        for lat, lon in zip(event_stns.latitude, event_stns.longitude)
    ]
    event_stns = GeoDataFrame(event_stns, geometry=geometry).set_crs(
        crs="EPSG:3857", allow_override=True
    )
    # adding geometry column
    event_stns_local = gpd.overlay(event_stns, target_counties, how="intersection")
    num_event_stns_local = len(event_stns_local)
    # subsetting for stations within county boundaries
    print(
        f"{num_event_stns_local} potential stations available for evaluation for {event_to_eval} event."
    )

    # Check if a specific_station is requested and return that one
    if specific_station is not None:
        eval_stns = event_stns[event_stns["era-id"] == specific_station]
        if len(eval_stns) == 0:
            raise ValueError(
                f"Station {specific_station} is not within the training/event dataset"
            )
        return eval_stns

    if subset != None:
        if num_event_stns_local <= subset:
            eval_stns = event_stns_local
        else:
            eval_stns = event_stns_local.sample(subset, replace=False)
            print(
                f"{subset} stations selected for evaluation for {event_to_eval} event!"
            )
    else:
        eval_stns = event_stns_local

    if return_stn_ids:
        print("Stations selected for evaluation:\n", list(eval_stns["era-id"]))

    return eval_stns


def id_all_flags(ds: xr.Dataset):
    """
    Prints all unique values of all eraqaqc flags

    Parameters
    ----------
    ds : xr.Dataset
        station data in xr format (not pd.DataFrame)

    Returns
    -------
    None
    """

    ds_vars = list(ds.keys())
    qc_vars = [i for i in ds_vars if "_eraqc" in i]
    if len(qc_vars) == 0:
        print(
            "Station has no eraqc variables -- please double check that this station has completed QA/QC!"
        )
    else:
        for var in qc_vars:
            print(var, np.unique(ds[var].data))

    return None


def event_info(
    event: str, alt_start_date: str | None = None, alt_end_date: str | None = None
) -> tuple[str, str]:
    """
    Utility function to return useful information for a specific designated event.

    Paramters
    ---------
    event : str
        name of event to evaluate
    alt_start_date : str
        date of different event, must be in format "YYYY-MM-DD"
    alt_end_date : str
        date of different event, must be in format "YYYY-MM-DD"

    Returns
    -------
    tuple[str, str]
        start and end dates for requested event

    To dos
    ------
    1. Alternative start / end date format check
    """

    start_date = {
        "santa_ana_wind": "2007-10-19",
        "winter_storm": "1990-12-20",
        "AR": "2017-01-16",
        "mudslide": "2018-01-05",
        "aug2020_heatwave": "2020-08-14",
        "sep2020_heatwave": "2020-09-05",
        "aug2022_heatwave": "2022-08-30",
        "offshore_wind": "2021-01-15",
        "alternative": alt_start_date,
    }

    end_date = {
        "santa_ana_wind": "2007-11-16",
        "winter_storm": "1990-12-24",
        "AR": "2017-01-20",
        "mudslide": "2018-01-09",
        "aug2020_heatwave": "2020-08-15",
        "sep2020_heatwave": "2020-09-08",
        "aug2022_heatwave": "2022-09-09",
        "offshore_wind": "2021-01-16",
        "alternative": alt_end_date,
    }

    event_start = start_date[event]
    event_end = end_date[event]

    return (event_start, event_end)


def event_subset(
    df: pd.DataFrame,
    event: str,
    buffer: int | None = 7,
    alt_start_date: str | None = None,
    alt_end_date: str | None = None,
) -> pd.DataFrame:
    """
    Subsets for the event itself + buffer around to identify event.

    Parameters
    ---------
    df: pd.DataFrame
        stationlist dataframe
    event : str
        name of event
    buffer : int, optional
        number of days to include as a buffer around event start/end date
    alt_start_date : str
        date of different event, must be in format "YYYY-MM-DD"
    alt_end_date : str
        date of different event, must be in format "YYYY-MM-DD"

    Returns
    -------
    event_sub : pd.DataFrame
        subset of stationlist within date range of event or alternative
    """

    print(
        f"Subsetting station record for event duration with {str(buffer)} day buffer..."
    )

    # set to searchable datetime
    df["time"] = pd.to_datetime(df["time"])
    # grab dates from lookup dictionary
    event_start, event_end = event_info(event, alt_start_date, alt_end_date)

    # subset for event dates + buffer
    datemask = (
        df["time"] >= (pd.Timestamp(event_start) - datetime.timedelta(days=buffer))
    ) & (df["time"] <= (pd.Timestamp(event_end) + datetime.timedelta(days=buffer)))
    event_sub = df.loc[datemask]

    return event_sub


def flags_during_event(subset_df: pd.DataFrame, var: str, event: str) -> list[str]:
    """
    Provides info on which flags were placed during event for evaluation

    Parameters
    ---------
    subset_df : pd.DataFrame

    var : str
        name of variable to assess flags
    event : str
        name of case study event

    Returns
    -------
    all_event_flags : list[str]
        all event flags set
    """

    event_flags = subset_df[var + "_eraqc"].unique()
    print(f"Flags set on {var} during {event} event: {event_flags}")
    all_event_flags = []
    for item in subset_df[var + "_eraqc"].unique():
        all_event_flags.append(item)

    return all_event_flags


def find_other_events(
    df, event_start, event_end, buffer=14, subset=None, return_stn_ids=True
):
    """
    Event finder not tied to specified case study events.

    Parameters
    ---------
    df : pd.DataFrame
        stationlist
    event_start : str
        start of event, format "YYYY-MM-DD"
    event_end : str
        end of event, format "YYYY-MM-DD"

    Returns
    -------
    eval_stns : pd.DataFrame
        subset of stations for other events of interest

    To dos
    ------
    1. Manual end date check no longer relevant, make sure stationlist passed is the correct updated version.
    2. Start / end date format check
    """

    print(
        f"Subsetting station record for event duration with {str(buffer)} day buffer..."
    )

    df["start-date"] = pd.to_datetime(df["start-date"])
    df["end-date"] = pd.to_datetime(df["end-date"])
    event_start = pd.to_datetime(event_start).tz_localize("UTC")
    event_end = pd.to_datetime(event_end).tz_localize("UTC")

    event_sub = df.loc[
        (df["start-date"] <= (event_start - datetime.timedelta(days=buffer)))
        & (df["end-date"] >= (event_end + datetime.timedelta(days=buffer)))
    ]

    # # exclude "manual check on end date" stations since we don't know when they actually end
    # event_sub = event_sub.loc[event_sub["notes"] != "manual check on end date"]

    # subset to make more manageable
    if subset != None:
        if len(event_sub) <= subset:
            eval_stns = event_sub
        else:
            eval_stns = event_sub.sample(subset, replace=False)
            print(f"{subset} stations selected for evaluation for comparison!")
    else:
        eval_stns = event_sub

    # return station ids for ease
    if return_stn_ids:
        print("Stations selected for evaluation:\n", list(eval_stns["era-id"]))

    return eval_stns


def latlon_to_mercator_cartopy(lat: float, lon: float) -> tuple[float, float]:
    """
    Converts lat/lon coordinates to mercator for plotting.

    Parameters
    ---------
    lat : float
        latitude
    lon : float
        longitude

    Returns
    -------
    tuple[float, float]
    """

    proj_latlon = CRS("EPSG:4326")
    proj_mercator = CRS("EPSG:3857")

    # Transform the coordinates
    transformer = Transformer.from_crs(proj_latlon, proj_mercator, always_xy=True)
    x, y = transformer.transform(lon, lat)

    return x, y


def stn_visualize(stn_id, stn_list, event_to_eval):
    """
    Produces simple map of station relevant to event boundary.

    Parameters
    ----------
    stn_id : str
        name of station
    stn_list : pd.DataFrame
        stationlist
    event_to_eval : str
        name of case study event

    Returns
    -------
    None
    """

    # grab station id info and reproject coords
    stn = stn_list.loc[stn_list["era-id"] == stn_id]
    lon, lat = stn.longitude.values[0], stn.latitude.values[0]
    x, y = latlon_to_mercator_cartopy(lat, lon)

    # figure set-up
    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.epsg(3857)})
    ax.coastlines()
    ax.add_feature(cf.BORDERS)
    ax.add_feature(cf.STATES, lw=0.5)

    ax.set_extent([lon + 1, lon - 1, lat - 1, lat + 1])

    # Obtain the limits of the plot
    x0, x1, y0, y1 = ax.get_extent()

    # Create a polygon with the limits of the plot
    polygon = Polygon(((x0, y0), (x0, y1), (x1, y1), (x1, y0)))

    # Use only the counties that overlap with the actual plot
    ca_county = gpd.read_file(CENSUS_SHP)
    counties = ca_county[ca_county.overlaps(polygon)]

    # Plot the counties' geometries
    for geometry in counties.geometry:
        ax.add_geometries(
            geometry.boundary,
            crs=ax.projection,
            facecolor="none",
            edgecolor="teal",
            lw=0.5,
        )

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    ax.plot(lon, lat, "ok", markersize=8, transform=ccrs.PlateCarree(), mfc="none")
    ax.plot(x, y, ".r", markersize=4)
    ax.annotate(f"{stn_id}", xy=(x, y), xytext=(x + 10, y + 10), fontsize=6)
    # station name
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(), draw_labels=["bottom", "left"], ls=":", lw=0.5
    )
    ax.set_title(f"{event_to_eval} evaluation \nat {stn_id}")

    return None


def event_plot(
    df: pd.DataFrame,
    var: str,
    event: str,
    alt_start_date: str | None = None,
    alt_end_date: str | None = None,
    dpi: int | None = None,
):
    """Produces timeseries of variables that have flags placed

    Parameters
    ----------
    df : pd.DataFrame
        stationlist
    var : str
        name of variable
    event : str
        name of case study event
    alt_start_date : str
        date of different event, must be in format "YYYY-MM-DD"
    alt_end_date : str
        date of different event, must be in format "YYYY-MM-DD"
    dpi : int
        figure resolution

    Returns
    -------
    None
    """

    fig, ax = plt.subplots(figsize=(10, 3))

    # plot all observations
    df.plot(
        ax=ax,
        x="time",
        y=var,
        marker=".",
        ms=4,
        lw=1,
        color="k",
        alpha=0.5,
        label="Cleaned data",
    )

    # plot event timeline
    event_start, event_end = event_info(event, alt_start_date, alt_end_date)
    ax.axvspan(event_start, event_end, color="red", alpha=0.1, label="{}".format(event))

    # ax.axhline(event_start, color='red', lw=2, alpha=0.25)
    # ax.axhline(event_end, color='red', lw=2, alpha=0.25)
    # ax.fill_between(x='time', 0, 1, where=y)

    # plot any flags placed by QA/QC
    if len(df[var + "_eraqc"].dropna().unique()) != 0:
        # identify flagged data, can handle multiple flags
        for flag in df[var + "_eraqc"].dropna().unique():
            flag_name = id_flag(flag)
            flag_str = 100 * len(df.loc[df[var + "_eraqc"] == flag, var]) / len(df)
            flag_label = f"{flag_str:.3f}% of data flagged by {flag_name}"

            flagged_data = df[~df[var + "_eraqc"].isna()]
            flagged_data.plot(
                x="time",
                y=var,
                ax=ax,
                marker="o",
                ms=7,
                lw=0,
                mfc="none",
                color="C3",
                label=flag_label,
            )

    legend = ax.legend(loc="upper left", prop={"size": 8})

    # plot aesthetics
    ylab, units, miny, maxy = _plot_format_helper(var)
    plt.ylabel(f"{ylab} [{units}]")
    plt.xlabel("")
    stn = df["station"].unique()[0]
    plt.title(
        f"QA/QC event evaluation: {event}: {stn}",
        fontsize=10,
    )

    return None
