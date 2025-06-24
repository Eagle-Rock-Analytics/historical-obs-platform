"""
figure_plots.py

Functions
---------
- get_hdp_colordict: Builds a dictionary of specified network colors for use in HDP figures.
- var_fullname: Returns the full name of variable.
- 

Intended Use
------------
"""

import boto3
import pandas as pd
import matplotlib.pyplot as plt
from io import BytesIO
import geopandas as gpd
import contextily as cx




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


## Plots independent of stage
## # of stations per network over time --> get_station_chart + new plot fn
## get_station_map 2 versions


## PULL

# Function: plot station chart
# Set update = True if you want to regenerate the primary station list, otherwise function pulls the existing file from AWS.
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

## CLEAN




def gdf_setup(var):
    # AWS set-up
    s3 = boto3.resource("s3")
    s3_cl = boto3.client("s3")
    bucket_name = "wecc-historical-wx"

    # Read in all stations
    obj = s3_cl.get_object(
        Bucket=bucket_name, Key="2_clean_wx/temp_clean_all_station_list.csv"
    )
    body = obj["Body"].read()
    df_all = pd.read_csv(BytesIO(body), encoding="utf8")

    # ------------------------------------------------------------------------------------------------------------
    # Make a geodataframe
    gdf = gpd.GeoDataFrame(
        df_all, geometry=gpd.points_from_xy(df_all.longitude, df_all.latitude)
    )
    gdf.set_crs(epsg=4326, inplace=True)  # Set CRS

    # Project data to match base tiles.
    gdf_wm = gdf.to_crs(epsg=3857)  # Web mercator

    # Read in geometry of WECC
    mar = gpd.read_file(
        "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"
    )
    ter = gpd.read_file(
        "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
    )
    wecc = gpd.GeoDataFrame(pd.concat([mar, ter]))

    # Use to clip stations
    wecc = wecc.to_crs(epsg=3857)
    gdf_wecc = gdf_wm.clip(wecc)
    gdf_wecc = gdf_wecc.sort_values(["network"])

    # Setting color
    gdf_wecc["Color"] = gdf_wecc["network"].map(color_dict)

    # Subsetting based on variable
    new_gdf = gdf_wecc.loc[gdf_wecc[str(var) + "_nobs"] > 0]

    return new_gdf


def single_var_map(var):
    # aws set-up
    s3 = boto3.resource("s3")
    s3_cl = boto3.client("s3")
    bucket_name = "wecc-historical-wx"

    # grab gdf per variable
    gdf = gdf_setup(var)

    # figure global settings
    fig, ax = plt.subplots(figsize=(9, 9))
    a = 1  # alpha
    ms = 2  # markersize

    # enforce WECC boundary, all plots regardless of variable selection
    ylim = [3.5e6, 8.5e6]  # lat
    xlim = [-1.53e7, -1.13e7]  # lon
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)

    # plot by network with correct legend color
    for ctype, data in gdf.groupby("network"):
        color = color_dict[ctype]
        data.plot(color=color, markersize=ms, alpha=a, ax=ax, label=ctype)

    cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)

    # legend
    l = ax.legend(loc="lower left", prop={"size": 8}, frameon=True)
    #     for i in range(len(l.legend_handles)):
    #         l.legend_handles[i]._sizes = [20] # sets size of points in legend to be more visible

    # set title
    vartitle = var_fullname(var)
    ax.set_title(vartitle, fontsize=10)
    ax.set_axis_off()

    # in notebook, show figure
    plt.show()

    # save to AWS
    img_data = BytesIO()
    plt.savefig(img_data, format="png", bbox_inches="tight")
    img_data.seek(0)

    bucket = s3.Bucket(bucket_name)
    bucket.put_object(
        Body=img_data,
        ContentType="image/png",
        Key="2_clean_wx/clean_station_{}.png".format(var),
    )

def combo_var_map(multi_var):
    # aws set up
    s3 = boto3.resource("s3")
    s3_cl = boto3.client("s3")
    bucket_name = "wecc-historical-wx"

    # figure global settings
    a = 1  # alpha
    ms = 2  # markersize
    n = 1

    # set up - determine how many vars are passed
    if len(multi_var) == 1:  # single var passed
        single_var_map(multi_var[0])

    elif len(multi_var) % 2 == 0:  # even number passed
        fig = plt.subplots(nrows=int((len(multi_var) / 2)), ncols=2, figsize=(16, 20))

        for var_to_plot in multi_var:
            ax = plt.subplot((len(multi_var) / 2), 2, n)

            # enforce WECC boundary, all plots regardless of variable selection
            ylim = [3.5e6, 8.5e6]  # lat
            xlim = [-1.53e7, -1.13e7]  # lon
            ax.set_ylim(ylim)
            ax.set_xlim(xlim)

            gdf = gdf_setup(var_to_plot)

            # plot by network with correct legend color
            for ctype, data in gdf.groupby("network"):
                color = color_dict[ctype]
                data.plot(color=color, markersize=ms, alpha=a, ax=ax, label=ctype)

            cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)

            # set title
            vartitle = var_fullname(var_to_plot)
            ax.set_title(vartitle, fontsize=12)
            ax.set_axis_off()

            # legend
            l = ax.legend(loc="lower left", prop={"size": 8}, frameon=True)
            #             for i in range(len(l.legend_handles)):
            #                 l.legend_handles[i]._sizes = [20] # sets size of points in legend to be more visible

            n = n + 1  # move to next axis

    elif len(multi_var) % 3 == 0:  # odd number passed
        fig = plt.subplots(nrows=int((len(multi_var) / 3)), ncols=3, figsize=(24, 20))

        for var_to_plot in multi_var:
            ax = plt.subplot(len(multi_var) / 3, 3, n)

            # enforce WECC boundary, all plots regardless of variable selection
            ylim = [3.5e6, 8.5e6]  # lat
            xlim = [-1.53e7, -1.13e7]  # lon
            ax.set_ylim(ylim)
            ax.set_xlim(xlim)

            gdf = gdf_setup(var_to_plot)

            # plot by network with correct legend color
            for ctype, data in gdf.groupby("network"):
                color = color_dict[ctype]
                data.plot(color=color, markersize=ms, alpha=a, ax=ax, label=ctype)

            cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)

            # set title
            vartitle = var_fullname(var_to_plot)
            ax.set_title(vartitle, fontsize=12)
            ax.set_axis_off()

            # legend
            l = ax.legend(loc="lower left", prop={"size": 8}, frameon=True)
            #             for i in range(len(l.legend_handles)):
            #                 l.legend_handles[i]._sizes = [20] # sets size of points in legend to be more visible

            # move to next axis
            n = n + 1

    # unifed single legend -- to do

    # fix white spacing
    plt.tight_layout()

    # Save to AWS
    img_data = BytesIO()
    plt.savefig(img_data, format="png", bbox_inches="tight")
    img_data.seek(0)

    bucket = s3.Bucket(bucket_name)
    bucket.put_object(
        Body=img_data,
        ContentType="image/png",
        Key="2_clean_wx/clean_station_{}.png".format("_".join(multi_var)),
    )

    return None


def plot_all_vars(var_list):
    # aws set up
    s3 = boto3.resource("s3")
    s3_cl = boto3.client("s3")
    bucket_name = "wecc-historical-wx"

    # subplot global settings
    a = 1  # alpha
    ms = 2  # markersize
    n = 1  # initialize plot counter

    fig = plt.subplots(nrows=3, ncols=5, figsize=(24, 18))

    for var_to_plot in var_list:
        ax = plt.subplot(3, 5, n)

        # enforce WECC boundary, all plots regardless of variable selection
        ylim = [3.5e6, 8.5e6]  # lat
        xlim = [-1.53e7, -1.13e7]  # lon
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

        gdf = gdf_setup(var_to_plot)

        # plot by network with correct legend color
        for ctype, data in gdf.groupby("network"):
            color = color_dict[ctype]
            data.plot(
                color=color, markersize=ms, alpha=a, ax=ax, label=ctype, legend=False
            )

        cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)

        # set title
        vartitle = var_fullname(var_to_plot)
        ax.set_title(vartitle, fontsize=12)
        ax.set_axis_off()

        # legend handling
        l = ax.legend(loc="lower left", prop={"size": 6}, frameon=True)
        #         for i in range(len(l.legend_handles)):
        #             l.legend_handles[i]._sizes = [1] # sets size of points in legend to be more visible

        # move ot next subplot
        n = n + 1

    # fix white spacing
    plt.tight_layout()

    plt.show()

    # Save to AWS
    img_data = BytesIO()
    plt.savefig(img_data, format="png", bbox_inches="tight")
    img_data.seek(0)

    bucket = s3.Bucket(bucket_name)
    bucket.put_object(
        Body=img_data,
        ContentType="image/png",
        Key="2_clean_wx/clean_station_allvars.png",
    )


## QAQC

## MERGE