"""
update_zarr_attributes.py

Script to add or update a specified attribute on all Zarr datasets referenced 
in an Intake-ESM catalog stored on S3, consolidating metadata to ensure the 
attribute is recognized by tools like Zarr and xarray.

Key features:
- Loads an Intake catalog from a public S3 URL.
- Iterates over each Zarr dataset path listed in the catalog.
- Opens each Zarr store in read/write mode.
- Adds or updates a customizable attribute key/value pair.
- Consolidates metadata for efficient cloud-based access.
- Uses multithreading with a progress bar to speed up processing.

Usage:
    python update_zarr_attributes.py

Configuration:
- Modify the `ATTR_KEY` and `ATTR_VALUE` constants to set the attribute to add or update.
- Adjust `max_workers` in `main()` to control parallelism depending on your network and system.

Requirements:
- Python packages: zarr, intake, tqdm
- Network access to the S3 bucket hosting the Zarr datasets

"""

import time
import zarr
import intake
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

ATTR_KEY = "DOI"
ATTR_VALUE = "https://zenodo.org/records/16370140"


def process_store(s3_path, attr_key, attr_value):
    """
    Open a Zarr store at the given path, add or update a specified attribute,
    and consolidate metadata to ensure the attribute is persisted and visible.

    Parameters
    ----------
    s3_path : str
        The path or URL to the Zarr store (e.g., S3 path).
    attr_key : str
        The attribute key to add or update in the Zarr store's root attributes.
    attr_value : Any
        The value to assign to the specified attribute key.

    Returns
    -------
    str
        The input `s3_path`, returned for optional logging or tracking.
    """
    # Open the Zarr store in read/write mode
    store = zarr.open(s3_path, mode="r+")

    # Add or update the custom attribute for the dataset
    store.attrs[attr_key] = attr_value

    # Consolidate metadata to update the .zmetadata file for cloud readers
    # This ensures that when you read the dataset in using xarray, the attribute shows up
    zarr.consolidate_metadata(store.store)
    return s3_path


def main():
    start_time = time.time()

    # Open Intake-ESM catalog from a public S3 URL
    cat = intake.open_esm_datastore(
        "https://cadcat.s3.amazonaws.com/histwxstns/era-hdp-collection.json"
    )
    paths_all = cat.df["path"].values

    # Number of worker threads to use (adjust as needed)
    max_workers = 10

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks to the thread pool
        futures = [
            executor.submit(process_store, path, ATTR_KEY, ATTR_VALUE)
            for path in paths_all
        ]

        # Use tqdm to show progress bar as futures complete
        for _ in tqdm(
            as_completed(futures), total=len(futures), desc="Updating Zarr attributes"
        ):
            pass

    elapsed = time.time() - start_time
    print(
        f"\nDone updating {len(paths_all)} datasets. Time elapsed: {elapsed:.2f} seconds."
    )


if __name__ == "__main__":
    main()
