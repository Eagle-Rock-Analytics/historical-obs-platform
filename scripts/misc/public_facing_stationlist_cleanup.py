"""
public_facing_stationlist_cleanup.py

Prepares a cleaned station list for public distribution by correcting elevation 
errors and creating standardized source IDs for merged weather station data.

Operations:
- Filters to successfully merged stations
- Corrects bad/missing ASOSAWOS elevation values using source data
- Creates source-id field (uses WBAN codes for ASOSAWOS, era-id suffix for others)
- Uses existing latitude and longitude columns to create a geometry column 
- Reads in Tiger US states shapefile and adds US state for each station 
- Exports subset of columns as CSV for public-facing applications
"""

import pandas as pd
import numpy as np
import geopandas as gpd

# s3 Paths
MERGE_LIST_PATH = "s3://wecc-historical-wx/4_merge_wx/all_network_stationlist_merge.csv"
ASOSAWOS_ISD_PATH = (
    "s3://wecc-historical-wx/1_raw_wx/ASOSAWOS/stationlist_ISD_ASOSAWOS.csv"
)

# CSV-related paths
CSV_OUTPUT_FILENAME = "historical_wx_stations.csv"
CSV_OUTPUT_FILEPATH = f"../../data/{CSV_OUTPUT_FILENAME}"


def main():
    # Read in tables from s3
    merge_df = pd.read_csv(MERGE_LIST_PATH, index_col=0)
    asosawos_df = pd.read_csv(ASOSAWOS_ISD_PATH, index_col=0)

    # Just get stations that completed merge step
    merge_df = merge_df[merge_df["merged"] == "Y"].reset_index(drop=True)

    # Add ASOSAWOS era_id
    asosawos_df["era-id"] = "ASOSAWOS_" + asosawos_df["ISD-ID"].str.replace("-", "")

    # Rename elevation column to match merge_df
    asosawos_df.rename(columns={"ELEV(M)": "elevation"}, inplace=True)

    # Check for ASOSAWOS bad or missing elevation in merge dataframe
    bad_elevation = [-30479.6952]
    asosawos_bad_elevation = merge_df[
        (merge_df["network"] == "ASOSAWOS")
        & ((merge_df["elevation"].isin(bad_elevation)) | (merge_df["elevation"].isna()))
    ]

    # Get correct elevation from ASOSAWOS raw csv by merging on era-id
    elevation_corrections = asosawos_df[["era-id", "elevation"]].copy()

    # Update elevation values in merge_df
    for idx, row in asosawos_bad_elevation.iterrows():
        era_id = row["era-id"]
        correct_elev = elevation_corrections[elevation_corrections["era-id"] == era_id][
            "elevation"
        ].values
        if len(correct_elev) > 0:
            merge_df.loc[idx, "elevation"] = correct_elev[0]

    # Get source ID by removing Network name and underscore
    merge_df["source-id"] = merge_df.apply(
        lambda row: row["era-id"].replace(row["network"] + "_", ""), axis=1
    )

    # Except, for Network=ASOSAWOS, the source_id should be ICAO column
    merge_df = merge_df.merge(
        asosawos_df[["era-id", "ICAO"]].astype({"ICAO": str}), on="era-id", how="left"
    )
    merge_df["source-id"] = np.where(
        merge_df["network"] == "ASOSAWOS",
        merge_df["ICAO"],
        merge_df["source-id"].astype(str),
    )
    merge_df.drop(columns=["ICAO"], inplace=True)

    # Add in geometry column so it can be used as a GeoDataFrame
    merge_df = gpd.GeoDataFrame(
        merge_df,
        geometry=gpd.points_from_xy(merge_df.longitude, merge_df.latitude),
        crs="EPSG:4326",
    )

    # Get US states shapefiles from US census Tiger dataset
    states = gpd.read_file(
        "https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_state_20m.zip"
    )
    states.to_crs(merge_df.crs, inplace=True)
    states.rename(columns={"STUSPS": "state"}, inplace=True)

    # Spatial join to find which state each station is in
    merge_df = gpd.sjoin(
        merge_df, states[["state", "geometry"]], how="left", predicate="within"
    )
    merge_df.drop(columns="index_right", inplace=True)

    # Convert geometry to WKT format
    # WKT (Well-Known Text) is a text format for representing vector geometries
    merge_df["geometry"] = merge_df.geometry.to_wkt()

    # Subset columns
    merge_df_public_facing = merge_df[
        [
            "era-id",
            "source-id",
            "network",
            "latitude",
            "longitude",
            "state",
            "elevation",
            "start-date",
            "end-date",
            "total_nobs",
            "geometry",
        ]
    ]

    # Export DataFrame to local csv
    merge_df_public_facing.to_csv(CSV_OUTPUT_FILEPATH, index=False)
    print(f"Output file to {CSV_OUTPUT_FILEPATH}")


if __name__ == "__main__":
    main()
