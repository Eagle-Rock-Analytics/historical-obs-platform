"""figure_utils.py

Resamples data to an hourly timestep according to standard conventions.

Functions
---------
- qaqc_flag_fcn: Used for resampling QAQC flag columns.
- _modify_infill: This function does two things: (1) Flags rows ('y' under 'standardized_infill') that were infilled by resampling in the hourly standardization
    process, where there were time gaps in the input dataframe. (2) Infills constant variables (ie those in "constant_vars") observations that were left empty because
    they were in a time gap.
- merge_hourly_standardization: Resamples meteorological variables to hourly timestep according to standard conventions.

Intended Use
------------
Script functions are used to standardize all variables to hourly temporal resolution as a part of the merge pipeline.
"""

from functools import reduce
import numpy as np
import pandas as pd
import logging
import inspect


# ---------------------------------------------------------
# get_station_list_paths(bucket_name, directory)


# stage 1
def get_station_list_paths(bucket_name, directory):
    # Set up variables
    s3 = boto3.resource("s3")
    s3_cl = boto3.client("s3")  # for lower-level processes
    get_last_modified = lambda obj: int(
        obj["LastModified"].strftime("%s")
    )  #  Write method to get last modified file

    # Read in all station lists.
    # Get list of folder prefixes
    response = s3_cl.list_objects_v2(
        Bucket=bucket_name, Prefix=directory, Delimiter="/"
    )

    networks = {"Network": [], "NetworkPath": [], "StationFile": []}

    for prefix in response["CommonPrefixes"]:  # For each folder path
        networkpath = prefix["Prefix"][:-1]
        networkname = networkpath.split("/")[-1]
        station_file = (
            s3.Bucket(bucket_name)
            .objects.filter(Prefix=networkpath + "/" + "stationlist_")
            .all()
        )
        if len(list(station_file)) == 1:  # If one item returned
            for item in station_file:  # Get station file
                networks["Network"].append(networkname)
                networks["NetworkPath"].append(networkpath)
                networks["StationFile"].append(item.key)
                break  # If more than one file of this format found in folder, just take the most recent.
        elif len(list(station_file)) == 0:  # If no items found using search above
            files = s3.Bucket(bucket_name).objects.filter(
                Prefix=networkpath + "/"
            )  # List all files in folder
            file = [
                file for file in files if "station" in file.key
            ]  # More general search for 'station'
            for (
                item
            ) in (
                file
            ):  # Keep all files found here. These files may be different (e.g. ISD ASOS/AWOS vs ASOS/AWOS station lists)
                networks["Network"].append(networkname)
                networks["NetworkPath"].append(networkpath)
                networks["StationFile"].append(item.key)
        elif (
            len(list(station_file)) > 1
        ):  # If more than one identically formatted station list returned (shouldn't happen), take most recent
            file = [
                obj.key
                for obj in sorted(
                    station_file, key=lambda x: x.last_modified, reverse=True
                )
            ]  # Sort station files by most recent edit
            networks["Network"].append(networkname)
            networks["NetworkPath"].append(networkpath)
            networks["StationFile"].append(
                file[0]
            )  # Add path to most recently changed file. (Tested and works)
        # Note method currently doesn't have a method for dealing with more than one normally formatted station file
    return networks


# stage 2
def get_station_list_paths(bucket_name, directory):
    # Set up variables
    s3 = boto3.resource("s3")
    s3_cl = boto3.client("s3")
    get_last_modified = lambda obj: int(
        obj["LastModified"].strftime("%s")
    )  #  Write method to get last modified file

    # Read in all station lists
    # Get list of folder prefixes
    response = s3_cl.list_objects_v2(
        Bucket=bucket_name, Prefix=directory, Delimiter="/"
    )

    networks = {"Network": [], "NetworkPath": [], "StationFile": []}

    for prefix in response["CommonPrefixes"]:  # For each folder path
        networkpath = prefix["Prefix"][:-1]
        networkname = networkpath.split("/")[-1]
        station_file = (
            s3.Bucket(bucket_name)
            .objects.filter(Prefix=networkpath + "/" + "stationlist_")
            .all()
        )

        if len(list(station_file)) == 1:  # If one item returned
            for item in station_file:  # Get station file
                networks["Network"].append(networkname)
                networks["NetworkPath"].append(networkpath)
                networks["StationFile"].append(item.key)
                break  # If more than one file of this format found in folder, just take the most recent

        elif len(list(station_file)) == 0:  # If no items found using search above
            files = s3.Bucket(bucket_name).objects.filter(
                Prefix=networkpath + "/"
            )  # List all files in folder
            file = [
                file for file in files if "station" in file.key
            ]  # More general search for 'station'

            # Keep all files found here. These files may be different (e.g. ISD ASOS/AWOS vs ASOS/AWOS station lists)
            for item in file:
                networks["Network"].append(networkname)
                networks["NetworkPath"].append(networkpath)
                networks["StationFile"].append(item.key)

        # If more than one identically formatted station list returned, take most recent
        elif len(list(station_file)) > 1:
            # sort station lists by last edit
            file_all = [obj.key for obj in station_file]

            # need to handle for subsetted CWOP stationfiles
            if networkname == "CWOP":
                networks["Network"].append(networkname)
                networks["NetworkPath"].append(networkpath)
                networks["StationFile"].append(
                    file_all[0]
                )  # Add path to first in list -- alphabetical

            else:
                file = [
                    obj.key
                    for obj in sorted(
                        station_file, key=lambda x: x.last_modified, reverse=True
                    )
                ]
                networks["Network"].append(networkname)
                networks["NetworkPath"].append(networkpath)
                networks["StationFile"].append(
                    file[0]
                )  # Add path to most recently changed file

    # Note: currently doesn't have a method for dealing with more than one normally formatted station file
    return networks


# ---------------------------------------------------------
# get_station_chart(bucket_name, directory)

# stage 1


def get_station_chart(bucket_name, directory, update=False):
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

    # Fix nas
    ## Filter out rows w/o start date
    ## Note here: we lose MARITIME and NDBC networks.
    # print(dffull[dffull['network']=="MARITIME"])
    subdf = dffull.loc[~dffull["start-date"].isnull()].copy()

    ## Filter out non-downloaded rows
    subdf = subdf.loc[subdf["pulled"] != "N"].copy()

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

    out = subdf.explode("period").pivot_table(
        values="name", index="network", columns="period", aggfunc="count", fill_value=0
    )
    # out.columns = out.columns.strftime('%b-%y')

    return out


# stage 2
def get_station_chart(bucket_name, directory, update=False):
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

    # Get period
    # Format dates in datetime format (this gets lost in import).
    dffull["start-date"] = pd.to_datetime(dffull["start-date"], utc=True)
    dffull["end-date"] = pd.to_datetime(dffull["end-date"], utc=True)

    # Fix nas
    # Filter out rows w/o start date
    # Note here: we lose MARITIME and NDBC networks
    # print(dffull[dffull['network']=="MARITIME"])
    subdf = dffull.loc[~dffull["start-date"].isnull()].copy()

    # Filter out non-cleaned rows
    subdf = subdf.loc[subdf["cleaned"] != "N"].copy()

    # Manually filter dates to >01-01-1980 and <today.
    # Timezones so far ignored here but we presume on the scale of month we can safely ignore them for the moment
    # Note!: This implicitly assumes stations w/o end date run until present
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

    out = subdf.explode("period").pivot_table(
        values="era-id",
        index="network",
        columns="period",
        aggfunc="count",
        fill_value=0,
    )
    # out.columns = out.columns.strftime('%b-%y')

    return out


# stage 3
def get_station_chart(bucket_name, directory):
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

    # Fix nas
    ## Filter out rows w/o start date
    ## Note here: we lose MARITIME and NDBC networks.
    # print(dffull[dffull['network']=="MARITIME"])
    subdf = dffull.loc[~dffull["start-date"].isnull()].copy()

    ## Filter out non-cleaned rows
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

    out = subdf.explode("period").pivot_table(
        values="era-id",
        index="network",
        columns="period",
        aggfunc="count",
        fill_value=0,
    )
    # out.columns = out.columns.strftime('%b-%y')

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
