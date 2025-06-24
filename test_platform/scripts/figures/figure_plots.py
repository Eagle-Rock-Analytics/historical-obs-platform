"""
figure_plots.py

Functions
---------
- get_hdp_colordict: Builds a dictionary of specified network colors for use in HDP figures.
- var_fullname: Returns the full name of variable.
- networks_over_time_barchart: 

Intended Use
------------
"""

import boto3
import pandas as pd
import matplotlib.pyplot as plt
from io import BytesIO
import geopandas as gpd
import contextily as cx

s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")

BUCKET_NAME = "wecc-historical-wx"
PULL_DIR = "1_raw_wx"
CLEAN_DIR = "2_clean_wx"
QAQC_DIR = "3_qaqc_wx"
MERGE_DIR = "4_merge_wx"
WECC_MAR = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"
WECC_TERR = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"


def _get_stage_dir(stage: str) -> str:
    """
    Based on option passed, returns the correct AWS folder. 
    
    Parameters
    ----------
    stage : str
        name of stage of development
    
    Returns
    -------
    stage_dir : str
        name of corresponding AWS folder
    """

    if stage == "pull":
        stage_dir = PULL_DIR
    if stage == "clean":
        stage_dir = CLEAN_DIR
    if stage == "qaqc":
        stage_dir = QAQC_DIR
    if stage == "merge":
        stage_dir = MERGE_DIR

    return stage_dir

def _get_stage_stnlist(stage: str) -> pd.DataFrame:
    """
    Retrieves the corresponding stage station list.
    
    Parameters
    ----------
    stage : str
        name of stage of development
        
    Returns
    -------
    stage_stnlist : pd.DataFrame
        corresponding stage station list    
    """

    # get corresponding dir
    stage_dir = _get_stage_dir(stage)

    # read in stationlist from corresponding stage
    stage_stnlist = pd.read_csv(f"s3://{BUCKET_NAME}/{stage_dir}/all_network_stationlist_{stage}.csv")

    return stage_stnlist



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

def thing():
    return None

def networks_over_time_barchart(df: pd.DataFrame, stage: str, save_fig: bool=False) -> matplotlib.fig:
    """
    Produces stacked bar chart for network coverage over time. 

    Parameters
    ----------
    df : pd.DataFrame
        stage specific dataframe to generate barchart
    stage : str
        stage of development to generate barchart. Options: pull, clean, qaqc, merge
    save_fig : bool
        saves figure to AWS

    Returns
    -------
    fig : matplotlib.fig
        figure object for display in nb
    """

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
        title=f"{stage} stations by network over time", # capitalize if possible
        ax=ax,
        x_compat=True,
        cmap="tab20c_r",
    )  
    # Get area plot
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))  # Fix legend
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

    # Annotate with number of stations
    num_stns = 0 #! 
    plt.annotate(
        f"Total # of {stage} stations: {num_stns}", 
        xy=(0.025, 0.95), 
        xycoords="axes fraction"
    )

    # Save to AWS
    if save_fig: 
        dir_to_save = _get_stage_dir(stage)
        img_data = BytesIO()
        plt.savefig(img_data, format="png")
        img_data.seek(0)

        s3 = boto3.resource("s3")
        bucket = s3.Bucket(BUCKET_NAME)
        bucket.put_object(
            Body=img_data, 
            ContentType="image/png", 
            Key=f"{dir_to_save}/{stage}_stations_over_time.png"
        )

    return fig


## PULL

# Function: plot station chart
def get_station_chart(stage, ):

    # retrieve correct stage stationlist
   stage_stnlist = _get_stage_stnlist(stage)

    # Format dates in datetime format (this gets lost in import).
    stage_stnlist["start-date"] = pd.to_datetime(stage_stnlist["start-date"], utc=True)
    stage_stnlist["end-date"] = pd.to_datetime(stage_stnlist["end-date"], utc=True)

    # Fix nas
    ## Filter out rows w/o start date
    ## Note here: we lose MARITIME and NDBC networks.
    # print(dffull[dffull['network']=="MARITIME"])
    subdf = stage_stnlist.loc[~stage_stnlist["start-date"].isnull()].copy()
   # TODO: make sure that end dates on SNOTEL / SCAN are now fixed

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




def gdf_setup(var, stage):
    """

    Parameters
    ----------
    var : str

    stage : str

    Returns
    -------
    gdf
    """

    # Read in all stations
    df_all = _get_stage_stnlist(stage)

    # Make a geodataframe
    gdf = gpd.GeoDataFrame(
        df_all, 
        geometry=gpd.points_from_xy(
            df_all.longitude, 
            df_all.latitude)
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