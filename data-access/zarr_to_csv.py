"""
zarr_to_csv.py

Convert a Zarr store from S3 into a local CSV file using xarray and pandas.

This script:
1. Loads a specified Zarr store from an S3 URI.
2. Squeezes the singleton 'station' dimension to avoid multi-index issues.
3. Converts the xarray Dataset to a pandas DataFrame.
4. Optionally extracts dataset-level metadata and writes it as a comment-style header.
5. Saves the DataFrame as a CSV file locally.

Usage
-----
1. Set the `NETWORK` and `STATION` variables near the top of the script to match
   valid entries in the data catalog.

2. (Optional) Update `EXPORT_DIR` to your desired local output directory.
   By default, the CSV will be saved to the current working directory.
   Ensure the export directory exists.

3. (Optional) Adjust `ADD_METADATA_HEADER` to include/exclude metadata in the CSV.

4. Run the script:

       python zarr_to_csv.py

5. To load the exported CSV in pandas, skipping metadata header lines, use:

       import pandas as pd
       csv_filepath = f"{EXPORT_DIR}{STATION}.csv"
       df = pd.read_csv(csv_filepath, comment='#', index_col=0)

Output
------
On success, a CSV file named `{STATION}.csv` will be saved to `EXPORT_DIR`,
including metadata (if enabled) and tabular data.
"""

import xarray as xr
from zarr.errors import GroupNotFoundError
from typing import List  # used for function type hints

# Network and station name
NETWORK = "MARITIME"
STATION = "MARITIME_AAMC1"

# Directory to export the CSV file
# Add a trailing slash if specifying a folder, e.g., EXPORT_DIR = 'data/'
EXPORT_DIR = ""  # Default: current working directory

# Include metadata as a header in the CSV?
# Set to False to export only the data without the metadata header
ADD_METADATA_HEADER = True


def _get_metadata_from_ds(ds: xr.Dataset, network: str, station: str) -> List[str]:
    """
    Extract metadata attributes from an xarray Dataset and return them as a list of
    comment-formatted header lines for CSV output.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset containing global attributes.
    network : str
        The station network name, used as a fallback if not present in `ds.attrs`.
    station : str
        The station identifier, used as a fallback if not present in `ds.attrs`.

    Returns
    -------
    List[str]
        A list of strings, each representing a metadata line prefixed with "#".
    """

    attrs_dict = {
        "raw_files_merged": "n/a",
        "anemometer_height_m": "n/a",
        "thermometer_height_m": "n/a",
        "barometer_elevation_m": "n/a",
        "comment": "n/a",
        "disclaimer": "n/a",
    }

    for attr in attrs_dict.keys():
        try:  # Get attribute from dataset
            attrs_dict[attr] = ds.attrs[attr]
        except KeyError:
            continue  # Keep default if not found

    # Define custom multi-line header (as comment lines)
    metadata_header = [
        f"# STATION: {station}",
        f"# NETWORK: {network}",
        f"# RAW FILES MERGED: {attrs_dict['raw_files_merged']}",
        f"# ANEMOMETER HEIGHT (METERS): {attrs_dict['anemometer_height_m']}",
        f"# THERMOMETER HEIGHT (METERS): {attrs_dict['thermometer_height_m']}",
        f"# BAROMETER ELEVATION (METERS): {attrs_dict['barometer_elevation_m']}",
        f"# COMMENT: {attrs_dict['comment']}",
        f"# DISCLAIMER: {attrs_dict['disclaimer']}",
    ]

    return metadata_header


def main():

    print("Starting script zarr_to_csv.py")
    print(
        f"Network: {NETWORK}\nStation: {STATION}\nExport directory: {EXPORT_DIR}\nAdd metadata header: {ADD_METADATA_HEADER}"
    )

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

    # Set csv filepath
    csv_filename = f"{STATION}.csv"  # Name of file. Include csv extention
    csv_full_filepath = f"{EXPORT_DIR}{csv_filename}"  # Full filepath

    # Add metadata header
    if ADD_METADATA_HEADER:
        # Retrieve and format metadata for csv header
        metadata_header = _get_metadata_from_ds(ds, NETWORK, STATION)

        # Write metadata header to file
        with open(csv_full_filepath, "w") as f:
            for line in metadata_header:
                f.write(f"{line}\n")

    # Export data to csv
    # Append to the existing csv if the metadata has already been written
    try:
        df.to_csv(
            csv_full_filepath,
            mode="a" if ADD_METADATA_HEADER else "w",
            index=True,
            header=True,
        )
        print(f"Data successfully exported to {csv_full_filepath}")
    except Exception as e:
        print(f"ERROR: Could not export data to {csv_full_filepath}")
        print(e)

    print("Script complete.")


if __name__ == "__main__":
    main()
