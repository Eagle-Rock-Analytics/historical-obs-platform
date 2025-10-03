"""
public_facing_stationlist_cleanup.py

Prepares a cleaned station list for public distribution by correcting elevation 
errors and creating standardized source IDs for merged weather station data.

Operations:
- Filters to successfully merged stations
- Corrects bad/missing ASOSAWOS elevation values using source data
- Creates source-id field (uses WBAN codes for ASOSAWOS, era-id suffix for others)
- Exports subset of columns as CSV for public-facing applications
"""

import pandas as pd
import numpy as np

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

    # Except, for Network=ASOSAWOS, the source_id should be WBAN column (airport id)
    merge_df = merge_df.merge(
        asosawos_df[["era-id", "WBAN"]].astype({"WBAN": str}), on="era-id", how="left"
    )
    merge_df["source-id"] = np.where(
        merge_df["network"] == "ASOSAWOS",
        merge_df["WBAN"],
        merge_df["source-id"].astype(str),
    )
    merge_df.drop(columns=["WBAN"], inplace=True)

    # Subset columns
    merge_df_public_facing = merge_df[
        [
            "era-id",
            "source-id",
            "latitude",
            "longitude",
            "elevation",
            "start-date",
            "end-date",
            "network",
            "total_nobs",
        ]
    ]

    # Export DataFrame to local csv
    merge_df_public_facing.to_csv(CSV_OUTPUT_FILEPATH, index=False)
    print(f"Output file to {CSV_OUTPUT_FILEPATH}")


if __name__ == "__main__":
    main()
