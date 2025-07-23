"""
update_zarr_attributes.py

Script to add or update a custom attribute on all Zarr datasets referenced in an Intake-ESM catalog,
and consolidate metadata to ensure the attribute is recognized by Zarr/xarray.

This script:
- Loads an Intake catalog from a public S3 URL.
- Iterates over each Zarr dataset path in the catalog.
- Opens each Zarr store in read/write mode.
- Adds or updates a specified attribute with a given value.
- Consolidates metadata for efficient cloud access.

Usage:
    python update_zarr_attributes.py
"""

import time
import zarr
import intake
from tqdm import tqdm


def main():
    start_time = time.time()
    print("Starting script update_zarr_attributes.py")

    # Open Intake-ESM catalog from a public S3 URL
    cat = intake.open_esm_datastore(
        "https://cadcat.s3.amazonaws.com/histwxstns/era-hdp-collection.json"
    )

    # Extract all Zarr store paths from the catalog dataframe
    paths_all = cat.df["path"].values

    # Attribute key and value to add/update
    attr_key = "DOI"  # change this to any attribute name you want
    attr_value = (
        "https://zenodo.org/records/16370140"  # change this to any value you want
    )

    # Loop over each Zarr path and add/update the attribute with progress bar
    print(f"Updating attributes for {len(paths_all)} zarr stores")
    for s3_path in tqdm(paths_all, desc="Updating Zarr attributes"):
        # Open the Zarr store in read/write mode
        store = zarr.open(s3_path, mode="r+")

        # Add or update the custom attribute for the dataset
        store.attrs[attr_key] = attr_value

        # Consolidate metadata to update the .zmetadata file for cloud readers
        # This ensures that when you read the dataset in using xarray, the attribute shows up
        zarr.consolidate_metadata(store.store)

    elapsed = time.time() - start_time
    print(
        f"\nDone updating {len(paths_all)} datasets. Time elapsed: {elapsed:.2f} seconds."
    )


if __name__ == "__main__":
    main()
