"""figure_utils.py



Functions
---------
- 
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

# Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes

# Set relative paths to other folders and objects in repository.
BUCKET_NAME = "wecc-historical-wx"
RAW_DIR = '1_raw_wx'
CLEAN_DIR = '2_clean_wx'
QAQC_DIR = "3_qaqc_wx"
MERGE_DIR = "4_merge_wx"
stations_csv_path = f"s3://{BUCKET_NAME}/{QAQC_DIR}/all_network_stationlist_qaqc.csv"

#! could create dictionary


# ---------------------------------------------------------
# get_station_chart(bucket_name, directory)

# final function
def get_station_chart(directory,stage):
    """
    Sums two input flag count dataframes. This is a helper function for sum_flag_counts().

    Parameters
    ----------
    stage: str
        "pull", "clean", "qaqc" or "merge"
    directory: string
        RAW_DIR, CLEAN_DIR, QAQC_FIR, or MERGE_DIR
    
    Returns
    -------
    summed_df: pd.DataFrame

    """

    ## Get station list
    station_list = pd.read_csv(f"s3://{BUCKET_NAME}/{directory}/all_network_stationlist_{stage}.csv")

    #! only in qaqc
    # read in qaqc training station list
    stns = pd.read_csv("../3_qaqc_data/qaqc_training_station_list.csv")

    station_list = station_list[station_list["era-id"].isin(stns["era-id"])]
    #! only in qaqc ^

    ## Get period

    # Format dates in datetime format (this gets lost in import).
    station_list["start-date"] = pd.to_datetime(station_list["start-date"], utc=True)
    station_list["end-date"] = pd.to_datetime(station_list["end-date"], utc=True)

    # Fix nas
    ## Filter out rows w/o start date
    # print(dffull[dffull['network']=="MARITIME"])
    subdf = station_list.loc[~station_list["start-date"].isnull()].copy()

    # ! which we filter out depends on the phase, add if here
    ## Filter out non-downloaded rows  #! raw only
    subdf = subdf.loc[subdf["pulled"] != "N"].copy()

    # Filter out non-cleaned rows #! clean and QAQC
    subdf = subdf.loc[subdf["cleaned"] != "N"].copy()

    # manually filter dates to >01-01-1980 and <today.
    # Timezones so far ignored here but we presume on the scale of month we can safely ignore them for the moment.
    # Note!: This implicitly assumes stations w/o end date run until present.
    subdf["start-date"] = subdf["start-date"].apply(
        lambda x: (
            x
            if x > datetime(1980, 1, 1, tzinfo=timezone.utc)
            else datetime(1980, 1, 1, tzinfo=timezone.utc)
        )
    )
    subdf["end-date"] = subdf["end-date"].apply(
        lambda x: (
            x
            if x < datetime.utcnow().replace(tzinfo=timezone.utc)
            else datetime.utcnow().replace(tzinfo=timezone.utc)
        )
    )

    # Get period of months for range of dates for each station
    subdf["period"] = [
        pd.period_range(*v, freq="M")
        for v in zip(subdf["start-date"], subdf["end-date"])
    ]

    subdf = subdf[subdf.period.str.len() > 0]
    subdf = subdf.reset_index(drop=True)

    ##! multiple options for this step

    #! 1. from raw phase function
    out = subdf.explode("period").pivot_table(
        values="name", index="network", columns="period", aggfunc="count", fill_value=0
    )

    #! 3. from clean phase function
    out = subdf.explode("period").pivot_table(
        values="era-id",
        index="network",
        columns="period",
        aggfunc="count",
        fill_value=0,
    )

    return out


# ---------------------------------------------------------
# get_station_map(bucket_name, directory, shapepath, update=False)


# stage 1
def get_station_map(bucket_name, directory, shapepath, update=False):
    s3 = boto3.resource("s3")
    s3_cl = boto3.client("s3")  # for lower-level processes
    if update == False:
        obj = s3_cl.get_object(
            Bucket=bucket_name, Key="1_raw_wx/temp_pull_all_station_list.csv"
        )
        body = obj["Body"].read()
        dffull = pd.read_csv(BytesIO(body), encoding="utf8")
    elif update == True:
        dffull = get_station_list(bucket_name, directory)

    # Get period

    # Format dates in datetime format (this gets lost in import).
    dffull["start-date"] = pd.to_datetime(dffull["start-date"], utc=True)
    dffull["end-date"] = pd.to_datetime(dffull["end-date"], utc=True)

    # Quality control.
    # Fix nas
    ## Filter out rows w/o start date
    #     subdf = dffull.loc[~dffull['start-date'].isnull()].copy()
    # Filter out rows without data between 1980 and now.
    #     subdf = subdf.loc[(subdf['start-date']<=datetime.utcnow().replace(tzinfo=timezone.utc)) & (subdf['end-date']>='1980-01-01')]

    subdf = dffull

    # Make a geodataframe.
    gdf = gpd.GeoDataFrame(
        subdf, geometry=gpd.points_from_xy(subdf.longitude, subdf.latitude)
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

    # Version 1 - full map
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

    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)
    bucket.put_object(
        Body=img_data, ContentType="image/png", Key="1_raw_wx/pull_station_map.png"
    )


# Run function - generate station map
def get_station_map(bucket_name, directory, shapepath, update=False):
    s3 = boto3.resource("s3")
    s3_cl = boto3.client("s3")
    if update == False:
        obj = s3_cl.get_object(
            Bucket=bucket_name, Key="2_clean_wx/temp_clean_all_station_list.csv"
        )
        body = obj["Body"].read()
        dffull = pd.read_csv(BytesIO(body), encoding="utf8")
    elif update == True:
        dffull = get_station_list(bucket_name, directory)

    # Format dates in datetime format (this gets lost in import).
    dffull["start-date"] = pd.to_datetime(dffull["start-date"], utc=True)
    dffull["end-date"] = pd.to_datetime(dffull["end-date"], utc=True)

    # ------------------------------------------------------------------------------------------------------------
    #     # Quality control (optional -- uncomment next 3 lines of code if desired)
    #     # Filter out rows w/o start date - this will remove NDBC, MARITIME, CW3E networks (no date coverage in stn list)
    #     subdf = dffull.loc[~dffull['start-date'].isnull()].copy()
    #     # Filter out rows without data between 1980 and now.
    #     subdf = subdf.loc[(subdf['start-date']<=datetime.utcnow().replace(tzinfo=timezone.utc))
    #                       & (subdf['end-date']>='1980-01-01')]
    subdf = dffull

    # ------------------------------------------------------------------------------------------------------------
    # Make a geodataframe.
    gdf = gpd.GeoDataFrame(
        subdf, geometry=gpd.points_from_xy(subdf.longitude, subdf.latitude)
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

    # ------------------------------------------------------------------------------------------------------------
    # Version 1 - full map
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

    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)
    bucket.put_object(
        Body=img_data, ContentType="image/png", Key="2_clean_wx/clean_station_map.png"
    )

    # ------------------------------------------------------------------------------------------------------------
    # Version 2 - only big networks
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

    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)
    bucket.put_object(
        Body=img_data,
        ContentType="image/png",
        Key="2_clean_wx/clean_station_map_min.png",
    )


# stage 3


# Run function - generate station map
def get_station_map(bucket_name, directory, shapepath):
    s3 = boto3.resource("s3")
    s3_cl = boto3.client("s3")  # for lower-level processes

    # read in cleaned station list
    obj = s3_cl.get_object(
        Bucket=bucket_name, Key="2_clean_wx/temp_clean_all_station_list.csv"
    )
    body = obj["Body"].read()
    dfall = pd.read_csv(BytesIO(body), encoding="utf8")

    # read in qaqc training station list
    stns = pd.read_csv("../3_qaqc_data/qaqc_training_station_list.csv")

    dffull = dfall[dfall["era-id"].isin(stns["era-id"])]

    # Get period

    # Format dates in datetime format (this gets lost in import).
    dffull["start-date"] = pd.to_datetime(dffull["start-date"], utc=True)
    dffull["end-date"] = pd.to_datetime(dffull["end-date"], utc=True)

    # Quality control.
    # Fix nas
    ## Filter out rows w/o start date
    subdf = dffull.loc[~dffull["start-date"].isnull()].copy()
    # Filter out rows without data between 1980 and now.
    subdf = subdf.loc[
        (subdf["start-date"] <= datetime.utcnow().replace(tzinfo=timezone.utc))
        & (subdf["end-date"] >= "1980-01-01")
    ]

    # Make a geodataframe.
    gdf = gpd.GeoDataFrame(
        subdf, geometry=gpd.points_from_xy(subdf.longitude, subdf.latitude)
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

    # Version 1 - full map
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

    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)
    bucket.put_object(
        Body=img_data,
        ContentType="image/png",
        Key="3_qaqc_wx/qaqc_training_station_map.png",
    )
