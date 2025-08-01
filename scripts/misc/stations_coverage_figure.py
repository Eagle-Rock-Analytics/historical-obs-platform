"""
stations_coverage_figure.py

Generate a map figure showing all finalized, standardized, and QAQC'd stations
from the Historic Data Platform (HDP). Stations are plotted with colors corresponding
to their source networks.

Functions
---------
- create_stations_gdf: Reads a CSV file containing station metadata and returns a GeoDataFrame
  with Point geometries for plotting.
- read_network_colors: Reads a network colormap file and returns a dictionary mapping network IDs
  to hex color codes.
- create_map: Generates and saves a station coverage map with stations colored by network and
  overlaid on a basemap.
- main: Main script workflow: checks output path, builds station GeoDataFrame, reads network colors,
  plots stations by network with a basemap, and saves the figure.
"""

import inspect
import os
from time import time

import contextily as ctx
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# from shapely.geometry import Point


NETWORK_COLORS_PATH = "../../data/network_colors.txt"
STATIONLIST_URI = "s3://wecc-historical-wx/4_merge_wx/all_network_stationlist_merge.csv"

FIGS_DIR = "../../figures/"  # Directory to save figures to (must already exist)
FIGURE_SAVE_PATH = f"{FIGS_DIR}stations_coverage_map.png"  # Path to save figure to


def create_stations_gdf(csv_filepath):
    """
    Create a GeoDataFrame from a CSV file containing station metadata.

    Reads a CSV file with latitude and longitude columns, constructs Point geometries,
    and returns a GeoDataFrame reprojected to Web Mercator for mapping.

    Parameters
    ----------
    csv_filepath : str
        Path to the CSV file containing station metadata, including 'latitude' and 'longitude' columns.

    Returns
    -------
    gdf : gpd.GeoDataFrame
        GeoDataFrame with station metadata and a geometry column in EPSG:3857.
    """

    print(f"{inspect.currentframe().f_code.co_name}: Starting...")

    # Read in table from AWS
    stations_df = pd.read_csv(csv_filepath, index_col=0)

    # Construct geometry from lat and lon
    stations_df["geometry"] = gpd.points_from_xy(
        stations_df["longitude"], stations_df["latitude"]
    )

    # Convert to GeoDataFrame and set coordinate reference system (CRS)
    gdf = gpd.GeoDataFrame(stations_df, geometry="geometry", crs="EPSG:4326")  # WGS84
    gdf = gdf[["network", "era-id", "geometry"]].rename(
        columns={"era-id": "station_id"}
    )

    # Reproject to Web Mercator for basemap
    gdf = gdf.to_crs(epsg=3857)

    print(f"{inspect.currentframe().f_code.co_name}: Completed successfully")

    return gdf


def read_network_colors(filepath):
    """
    Read a network colormap file and return a dictionary mapping network IDs to hex colors.

    Parameters
    ----------
    filepath : str
        Path to the text file containing network colors. Each line should have two values:
        network and hex color (without the leading '#') separated by whitespace.

    Returns
    -------
    dict
        Dictionary where keys are network IDs (str) and values are hex color strings (e.g., '#d9d9d9').
    """
    print(f"{inspect.currentframe().f_code.co_name}: Starting...")
    color_dict = {}
    with open(filepath, "r") as f:
        for line in f:
            key, val = line.split()
            color_dict[key] = f"#{val}"

    print(f"{inspect.currentframe().f_code.co_name}: Completed successfully")
    return color_dict


def create_map(stations_gdf, colors_dict, fig_filepath, dpi=300):
    """
    Generate and save a station coverage map with stations colored by network.

    Parameters
    ----------
    stations_gdf : geopandas.GeoDataFrame
        GeoDataFrame containing station data with a 'network' column and Point geometries
        projected to Web Mercator (EPSG:3857).
    colors_dict : dict
        Dictionary mapping network IDs (str) to hex color codes (str), e.g. {'ASOSAWOS': '#d9d9d9'}.
    fig_filepath : str
        File path (including filename) where the resulting figure PNG will be saved.
    dpi : int, optional
        Resolution (dots per inch) for the saved figure. Default is 300.

    Returns
    -------
    None
        Saves the plot image to the specified path and prints progress messages.
    """

    print(f"{inspect.currentframe().f_code.co_name}: Starting...")

    # Set up figure
    fig, ax = plt.subplots(figsize=(10, 12))

    for network in stations_gdf["network"].unique():
        subset = stations_gdf[stations_gdf["network"] == network]
        subset.plot(
            ax=ax,
            color=colors_dict.get(network, "#000000"),  # Default color is black
            markersize=8,
            alpha=0.8,
            legend=False,  # disable legend in each plot call
        )

    # Create legend handles manually
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label=net,
            markerfacecolor=colors_dict[net],
            markersize=8,
        )
        for net in stations_gdf["network"].unique()
    ]

    ax.legend(
        handles=legend_handles,
        title="Source Network",
        loc="upper right",
        ncol=3,
        fontsize=11,
        title_fontsize=13,
        frameon=True,
    )

    # Add basemap
    ctx.add_basemap(ax, source=ctx.providers.CartoDB.Positron)

    # Set title and remove axis
    ax.set_title("Historic Data Platform Station Coverage", fontsize=16, pad=10)
    ax.axis("off")

    # Save
    plt.tight_layout()
    plt.savefig(fig_filepath, dpi=dpi)
    print(f"Figure saved to: {fig_filepath}")

    print(f"{inspect.currentframe().f_code.co_name}: Completed successfully")


def main():

    print("Starting script stations_coverage_figure.py")

    # Check if figure directory exists
    if not os.path.isdir(FIGS_DIR):
        raise NotADirectoryError(f"Figure directory does not exist: {FIGS_DIR}")

    # Create GeoDataFrame with station_id, network, and geometry
    stations_gdf = create_stations_gdf(STATIONLIST_URI)

    # Get colors corresponding to each network
    colors_dict = read_network_colors(NETWORK_COLORS_PATH)

    # Generate figure
    create_map(stations_gdf, colors_dict, FIGURE_SAVE_PATH, dpi=300)

    print("Script complete.")


if __name__ == "__main__":
    main()
