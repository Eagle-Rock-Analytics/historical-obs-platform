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

# CA county shapefile
CENSUS_SHP = "s3://wecc-historical-wx/0_maps/ca_counties/CA_Counties.shp"

# Event areas, v1 case studies -- fix dates
MUDSLIDE_1998_COUNTIES = ["Santa Barbara"]
HEATWAVE_2022_COUNTIES = [
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
SANTAANA_WIND_2007_COUNTIES = [
    "Los Angeles",
    "Orange",
    "San Diego",
    "San Bernardino",
    "Riverside",
    "Ventura",
    "Santa Barbara",
]
WINTERSTORM_1980_COUNTIES = [
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


## Event station finders
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


def _add_df_geometry(df, event):
    """ """

    # Convert df to gdf
    stns_gdf = gpd.GeoDataFrame(
        df, geometry=gpd.points_from_xy(df.longitude, df.latitude, crs="EPSG:4326")
    )

    # Convert to station CRS
    ca_counties = gpd.read_file(CENSUS_SHP)
    ca_counties = ca_counties.to_crs(stns_gdf.crs)

    # Build event area by counties
    if event == "aug2022_heatwave":
        event_counties = HEATWAVE_2022_COUNTIES
    elif event == "santa_ana_wind":
        event_counties = SANTAANA_WIND_2007_COUNTIES
    elif event == "mudslide":
        event_counties = MUDSLIDE_1998_COUNTIES
    elif event == "winter_storm":
        event_counties = WINTERSTORM_1980_COUNTIES

    # Define event geometry
    event_geom = ca_counties[ca_counties["NAME"].isin(event_counties)]

    # Filter down to stations that are in the event area
    stns_gdf["intersects"] = stns_gdf.intersects(event_geom.unary_union)
    # See which stations intersect with the event polygon

    # Get just those stations, drop the others
    event_stns = stns_gdf[stns_gdf["intersects"] == True].reset_index(drop=True)

    return event_stns


def find_event_stations(
    event: str,
    event_start: str,
    event_end: str,
    buffer: int = 14,
    subset: int | None = None,
    return_stn_ids: bool = True,
):
    """
    Event finder not tied to specified case study events.

    Parameters
    ---------
    event : str
        name of event, used to geographically subset stations
    event_start : str
        start of event, format "YYYY-MM-DD"
    event_end : str
        end of event, format "YYYY-MM-DD"
    buffer : int
        number of days to include around event start / end date
    subset : bool, optional
        number of stations to randomly subset for to manage size
    return_stn_ids : bool, optional
        option to return a list of stations with coverage during designated event


    Returns
    -------
    eval_stns : pd.DataFrame
        subset of stations for other events of interest
    """

    # Read in all stations list
    df = pd.read_csv(
        "s3://wecc-historical-wx/4_merge_wx/all_network_stationlist_merge.csv"
    )

    # Narrow list to stns that completed QAQC and standardization pipeline
    df = df.loc[df["merged"] == "Y"]

    print(
        f"Subsetting station record for event duration with {str(buffer)} day buffer..."
    )

    # Narrow list to stns that have temporal coverage during event
    df["start-date"] = pd.to_datetime(df["start-date"])
    df["end-date"] = pd.to_datetime(df["end-date"])
    event_start = pd.to_datetime(event_start).tz_localize("UTC")
    event_end = pd.to_datetime(event_end).tz_localize("UTC")

    event_sub = df.loc[
        (df["start-date"] <= (event_start - datetime.timedelta(days=buffer)))
        & (df["end-date"] >= (event_end + datetime.timedelta(days=buffer)))
    ]

    print(f"Number of stations with temporal coverage during event: {len(event_sub)}.")

    # Narrow list to stns that have spatial coverage during event
    event_sub = _add_df_geometry(event_sub, event)
    print(
        f"Number of stations with spatial coverage during event: {len(event_sub)}. \n"
    )

    # Subset to make more manageable
    if subset != None:
        if len(event_sub) <= subset:
            eval_stns = event_sub
        else:
            eval_stns = event_sub.sample(subset, replace=False)
            print(f"{subset} stations selected for evaluation for comparison!")
    else:
        eval_stns = event_sub

    # Return station ids for ease
    if return_stn_ids:
        print("Stations selected for evaluation:\n", list(eval_stns["era-id"]))

    return eval_stns


def find_event_vars(df: pd.DataFrame, var: list[str]):
    """ """

    # Subset based on inputted variables
    # Should be able to do multiple variables

    var_map = [v + "_nobs" for v in var]

    # One variable
    if len(var) == 1:
        df_var = df.loc[df[var_map[0]] > 0]

    else:
        print("TBD")

    return df_var


## Support
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


## Visualizations
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
