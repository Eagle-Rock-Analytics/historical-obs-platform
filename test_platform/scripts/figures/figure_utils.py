"""figure_utils.py

This script contains functions for visualizing stations over time in a map and chart.

Functions
---------
- get_station_chart(): prepares station list data to generate a chart of stations acquired over time, for a given phase
- plot_chart(): plots the output from get_station_chart(), exporting it to AWS
- get_station_map_v1(): generates and exports a map of station locations from the station list of the input phase
- get_station_map_v2(): generates and exports a map of large network station locations from the station list of the input phase

Intended Use
------------
Script functions used in visualization notebooks to generate figures
"""

import boto3
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import matplotlib
from io import BytesIO, StringIO
from datetime import datetime, timezone, date
import geopandas as gpd
import contextily as cx
from collections import defaultdict

# Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes

# Set relative paths to other folders and objects in repository.
BUCKET_NAME = "wecc-historical-wx"
RAW_DIR = "1_raw_wx"
CLEAN_DIR = "2_clean_wx"
QAQC_DIR = "3_qaqc_wx"
MERGE_DIR = "4_merge_wx"
phase_dict = {"pull": RAW_DIR, "clean": CLEAN_DIR, "qaqc": QAQC_DIR, "merge": MERGE_DIR}


def get_hdp_colordict() -> dict:
    """
    Builds a dictionary of specified network colors for use in HDP figures,
    using "network_colors.txt", which is a customized combo of the
    "tab20c_r" and "tab20b" matplotlib colormaps.

    Parameters
    ----------
    None

    Returns
    -------
    color_dict : dict
        dictionary of network names to designated colors
    """

    # initialize color dictionary
    color_dict = {}

    # read through network_colors text file and assign each network to its designated color
    with open("network_colors.txt") as f:
        for line in f:
            (key, val) = line.split()
            color_dict[key] = str("#") + str(val)

    return color_dict


def var_fullname(var: str) -> str:
    """
    Returns the full name of variable.

    Parameters
    ----------
    var : str
        name of variable

    Returns
    --------
    var_title : str
        long name of variable
    """

    if var == "tas":
        var_title = f"Air temperature ({var})"
    if "tdps" in var:
        var_title = f"Dewpoint temperature ({var})"
    elif var == "hurs":
        var_title = f"Relative humidity ({var})"
    elif var == "rsds":
        var_title = f"Radiation ({var})"
    elif var == "sfcwind":
        var_title = f"Surface wind speed ({var})"
    elif var == "sfcwind_dir":
        var_title = f"Surface wind direction ({var})"
    elif "pr" in var:
        var_title = f"Precipitation ({var})"
    elif "ps" in var and "td" not in var:
        var_title = f"Air pressure ({var})"

    return var_title


def get_station_chart(phase: str) -> tuple[pd.DataFrame, int] | None:
    """
    Prepares station list data to generate a chart of stations acquired over time, for a given phase.

    Parameters
    ----------
    phase: str
        the pipeline phase - "pull", "clean", "qaqc" or "merge"

    Returns
    -------
    out: pd.DataFrame
        station list modified for charting

    """

    if phase not in ["pull", "clean", "qaqc", "merge"]:
        print(f"invalid phase:{phase}")
        return None

    # Get station list
    directory = phase_dict[phase]
    station_list = pd.read_csv(
        f"s3://{BUCKET_NAME}/{directory}/all_network_stationlist_{phase}.csv"
    )

    # Get period
    # Format dates in datetime format (this gets lost in import).
    station_list["start-date"] = pd.to_datetime(station_list["start-date"], utc=True)
    station_list["end-date"] = pd.to_datetime(station_list["end-date"], utc=True)

    # Manually filter dates to >01-01-1980 and <today.
    # Timezones so far ignored here but we presume on the scale of month we can safely ignore them for the moment.
    # Note!: This implicitly assumes stations w/o end date run until present.
    station_list["start-date"] = station_list["start-date"].apply(
        lambda x: (
            x
            if x > datetime(1980, 1, 1, tzinfo=timezone.utc)
            else datetime(1980, 1, 1, tzinfo=timezone.utc)
        )
    )
    station_list["end-date"] = station_list["end-date"].apply(
        lambda x: (
            x
            if x < datetime.utcnow().replace(tzinfo=timezone.utc)
            else datetime.utcnow().replace(tzinfo=timezone.utc)
        )
    )

    # Get period of months for range of dates for each station
    station_list["period"] = [
        pd.period_range(*v, freq="M")
        for v in zip(station_list["start-date"], station_list["end-date"])
    ]

    subdf = station_list[station_list.period.str.len() > 0]
    subdf = subdf.reset_index(drop=True)

    # identify number of stations in each phase
    if phase == "qaqc":
        phased = "qaqc"
    elif phase == "merge":
        phased = "merged"
    else:
        phased = phase + "ed"
    phase_stn_count = len(subdf.loc[subdf[phased] != "N"])

    if phase == "pull":
        out = subdf.explode("period").pivot_table(
            values="name",
            index="network",
            columns="period",
            aggfunc="count",
            fill_value=0,
        )
    else:
        out = subdf.explode("period").pivot_table(
            values="era-id",
            index="network",
            columns="period",
            aggfunc="count",
            fill_value=0,
        )

    return out, phase_stn_count


def plot_chart(phase: str) -> None:
    """
    Plots the output from get_station_chart(), exporting it to AWS.

    Parameters
    ----------
    phase: str
        the pipeline phase - "pull", "clean", "qaqc" or "merge"

    Returns
    -------
    None
    """

    out, phase_stn_count = get_station_chart(phase)
    # Plot
    outt = out.T.reset_index()

    # Fix time component
    outt["date"] = outt["period"].astype(str)
    outt["date"] = pd.to_datetime(outt["date"])

    # Plot parameters
    plt.rcParams["figure.figsize"] = [7.50, 3.50]
    plt.rcParams["figure.autolayout"] = True
    plt.rcParams["figure.facecolor"] = "white"

    # Subplot parameters
    fig, ax = plt.subplots(figsize=(8, 6))
    outt.plot.area(
        x="date",
        title="Stations by network over time",
        ax=ax,
        x_compat=True,
        cmap="tab20c_r",
    )  # Get area plot
    ax.legend(loc="upper left", bbox_to_anchor=(0.01, 0.95), ncol=2)
    ax.tick_params(labelcolor="black", labelsize="medium", width=3)
    ax.set_facecolor("w")
    ax.set_xlabel("Date")
    ax.set_ylabel("Number of stations")

    # Change axis bounds
    ax.set_xlim([date(1980, 1, 1), date(2022, 8, 1)])

    # Change tick marks
    ax.minorticks_on()
    ax.xaxis.set_major_locator(matplotlib.dates.YearLocator(3))
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
    ax.xaxis.set_minor_locator(matplotlib.dates.YearLocator(1))

    # Change y ticks
    plt.locator_params(axis="y", nbins=12)
    ax.yaxis.get_ticklocs(minor=True)

    # Set x axis labels
    plt.subplots_adjust(left=0.2, bottom=0.2, top=0.8, right=0.8)

    # Annotate text for total number
    plt.annotate(
        f"Total # of {phase} stations: {phase_stn_count}",
        xy=(0.025, 0.95),
        xycoords="axes fraction",
    )

    # Save to AWS
    img_data = BytesIO()
    plt.savefig(img_data, format="png")
    img_data.seek(0)

    bucket = s3.Bucket(BUCKET_NAME)
    export_folder = phase_dict[phase]
    export_key = f"{export_folder}/{phase}_stations_over_time.png"
    bucket.put_object(
        Body=img_data,
        ContentType="image/png",
        Key=export_key,
    )

    return None


def get_station_map_v1(phase: str, shapepath: str) -> None:
    """
    Generates and exports a map of station locations from the station list of the input phase.

    Parameters
    ----------
    phase: str
        the pipeline phase - "pull", "clean", "qaqc" or "merge"
    shapepath: string
        path to shapefile ofr 2021 US state borders

    Returns
    -------
    None

    """
    if phase not in ["pull", "clean", "qaqc", "merge"]:
        print(f"invalid phase:{phase}")
        return None

    ## Get station list
    directory = phase_dict[phase]
    station_list = pd.read_csv(
        f"s3://{BUCKET_NAME}/{directory}/all_network_stationlist_{phase}.csv"
    )

    ## Get period

    # Format dates in datetime format (this gets lost in import).
    station_list["start-date"] = pd.to_datetime(station_list["start-date"], utc=True)
    station_list["end-date"] = pd.to_datetime(station_list["end-date"], utc=True)

    # Make a geodataframe.
    gdf = gpd.GeoDataFrame(
        station_list,
        geometry=gpd.points_from_xy(station_list.longitude, station_list.latitude),
    )
    gdf.set_crs(epsg=4326, inplace=True)  # Set CRS

    # Project data to match base tiles.
    gdf_wm = gdf.to_crs(epsg=3857)  # Web mercator

    # Read in geometry of continental US.
    us = gpd.read_file(shapepath)

    # Remove territories, AK, HI
    rem_list = ["HI", "AK", "MP", "GU", "AS", "PR", "VI"]
    us = us.loc[us.STUSPS.isin(rem_list) == False]

    # Use to clip stations
    us = us.to_crs(epsg=3857)
    gdf_us = gdf_wm.clip(us)

    # Plot
    ax = gdf_us.plot(
        "network",
        figsize=(15, 15),
        alpha=1,
        markersize=3,
        legend=True,
        cmap="nipy_spectral",
    )
    cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)
    ax.set_axis_off()

    # Save to AWS
    img_data = BytesIO()
    plt.savefig(img_data, format="png")
    img_data.seek(0)

    bucket = s3.Bucket(BUCKET_NAME)
    export_folder = phase_dict[phase]
    export_key = f"{export_folder}/{phase}_station_map.png"
    bucket.put_object(
        Body=img_data,
        ContentType="image/png",
        Key=export_key,
    )

    return None


def get_station_map_v2(phase: str, shapepath: str) -> None:
    """
    Generates and exports a map of large network station locations from the station list of the input phase.

    Rules
    -----
    1. Map networks with >100 stations, binning all other station under "misc"

    Parameters
    ----------
    phase: str
        the pipeline phase - "pull", "clean", "qaqc" or "merge"
    shapepath: string
        path to shapefile ofr 2021 US state borders

    Returns
    -------
    None

    """
    if phase not in ["pull", "clean", "qaqc", "merge"]:
        print(f"invalid phase:{phase}")
        return None

    ## Get station list
    directory = phase_dict[phase]
    station_list = pd.read_csv(
        f"s3://{BUCKET_NAME}/{directory}/all_network_stationlist_{phase}.csv"
    )

    ## Get period

    # Format dates in datetime format (this gets lost in import).
    station_list["start-date"] = pd.to_datetime(station_list["start-date"], utc=True)
    station_list["end-date"] = pd.to_datetime(station_list["end-date"], utc=True)

    # Make a geodataframe.
    gdf = gpd.GeoDataFrame(
        station_list,
        geometry=gpd.points_from_xy(station_list.longitude, station_list.latitude),
    )
    gdf.set_crs(epsg=4326, inplace=True)  # Set CRS

    # Project data to match base tiles.
    gdf_wm = gdf.to_crs(epsg=3857)  # Web mercator

    # Read in geometry of continental US.
    us = gpd.read_file(shapepath)

    # Remove territories, AK, HI
    rem_list = ["HI", "AK", "MP", "GU", "AS", "PR", "VI"]
    us = us.loc[us.STUSPS.isin(rem_list) == False]

    # Use to clip stations
    us = us.to_crs(epsg=3857)
    gdf_us = gdf_wm.clip(us)

    # Sort stations by number of networks
    gdf_us["network_count"] = gdf_us.groupby("network")["network"].transform(
        "count"
    )  # Add network count column.

    # If <100 stations, change to "misc"
    gdf_us.loc[gdf_us["network_count"] < 100, "network"] = "Misc"

    # Plot
    ax = gdf_us.plot(
        "network",
        figsize=(15, 15),
        alpha=1,
        markersize=3,
        legend=True,
        cmap="nipy_spectral",
    )
    cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)
    ax.set_axis_off()

    # Save to AWS
    img_data = BytesIO()
    plt.savefig(img_data, format="png")
    img_data.seek(0)

    bucket = s3.Bucket(BUCKET_NAME)
    export_folder = phase_dict[phase]
    export_key = f"{export_folder}/{phase}_station_map_min.png"
    bucket.put_object(
        Body=img_data,
        ContentType="image/png",
        Key=export_key,
    )

    return None
