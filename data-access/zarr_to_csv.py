"""
zarr_to_csv.py

Convert a Zarr store in S3 to a local CSV file using xarray.

This script:
1. Loads a specified Zarr store from an S3 URI.
2. Squeezes singleton 'station' dimension to avoid multi-index issues.
3. Converts the xarray dataset to a pandas DataFrame.
4. Exports the DataFrame as a CSV to a local directory.
"""

import xarray as xr
import pandas as pd
from zarr.errors import GroupNotFoundError

# Network and station name
NETWORK = "MARITIME"
STATION = "MARITIME_AAMC1"

# Where to export the csv file
# Add a trailing slash if EXPORT_DIR is a folder; i.e. EXPORT_DIR = 'data/'
EXPORT_DIR = ""  # Export to current directory


def main():

    print("Starting script zarr_to_csv.py")

    # Read in zarr store from s3 bucket
    s3_uri = f"s3://cadcat/histwxstns/{NETWORK}/{STATION}.zarr/"  # Build S3 path
    try:
        ds = xr.open_zarr(s3_uri)
    except GroupNotFoundError:
        print(f"ERROR: Zarr store not found in S3 at path: {s3_uri}")
        return None

    # Convert xr.Dataset object to pd.DataFrame
    ds_squeezed = ds.squeeze()
    df = ds_squeezed.to_dataframe()

    # Export to CSV
    csv_filename = f"{STATION}.csv"  # Name of file. Include csv extention
    csv_full_filepath = f"{EXPORT_DIR}{csv_filename}"  # Full filepath
    try:
        df.to_csv(csv_full_filepath, index=True)
        print(f"Data successfully exported to {csv_full_filepath}")
    except Exception as e:
        print(f"ERROR: Could not export data to {csv_full_filepath}")
        print(e)

    print("Script complete.")


if __name__ == "__main__":
    main()
