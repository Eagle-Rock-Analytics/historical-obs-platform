""" add_stations_to_csv.py 

Add station names to temp_clean_all_station_list.csv 
The script reads all the zarr files in s3 bucket 2 (wecc-historical-wx/2_clean_wx/STATION_NAME/) 
Then, retrieves the filename and elevation from each zarr and appends it (writes it) to the existing csv
Lastly, the updated csv is outputted to s3 or locally

NOTE: this script only works with zarr files (at the moment)

Author: Nicole Keeney 
Creation Date: Nov 27, 2024
Modification History: n/a 

"""

## ================= PACKAGE IMPORTS =================
import numpy as np
import pandas as pd
import xarray as xr
import warnings  # silence warnings in function call
from qaqc_utils import get_filenames_in_s3_folder, progressbar

## ================= GLOBAL VARIABLES =================
# s3 info
bucket = "wecc-historical-wx"
folder = "2_clean_wx"
network = "VALLEYWATER"

# csv info
csv_filename = "temp_clean_all_station_list.csv"
csv_local_path = csv_filename
csv_s3_filepath = "s3://{0}/{1}".format(bucket, csv_filename)


## ================= HELPER FUNCTION(S) =================
def _warn_if_more_than_one_value_in_array(l, name=""):
    """Print warning if there is more than one item in the array

    Parameters
    ----------
    l: np.array
    name: str, optional
        Name of variable to print to console (i.e. "longitude")

    Returns
    ------
    None

    """
    if len(l) > 1:
        print(
            "WARNING: More than one value found for {}. Using the first value from the following array: ".format(
                name
            )
        )
        print(repr(l))
    else:
        return None


## ================= MAIN FUNCTIONS =================
def add_stations_to_csv():
    print("csv will be saved to s3 path: {0}".format(csv_s3_filepath))

    # Read in file
    stations_df = pd.read_csv(csv_local_path)

    # Sort into alphabetical order
    stations_df = stations_df.sort_values("era-id")

    # Get allllll the files in the folder (including nested folders-- zarr stores are "folders" of many files as well)
    network_filenames = get_filenames_in_s3_folder(
        bucket=bucket, folder="{0}/{1}".format(folder, network)
    )

    # We just want to get the top directory for each station, i.e. "VALLEYWATER_6001.zarr/"
    network_filenames = [
        file.split(".zarr/")[0] + ".zarr/"
        for file in network_filenames
        if ".zarr/" in file
    ]

    # Since each station has a bunch of individual zarr stores, the split() function returns many copies of the same string
    # We just want one filename per station
    network_filenames = [
        x for i, x in enumerate(network_filenames) if x not in network_filenames[:i]
    ]
    if len(network_filenames) < 1:
        raise ValueError(
            "No zarrs found at path s3://{0}/{1}/{2}/\nbucket: {0}\nfolder: {1}\nnetwork: {2}\n".format(
                bucket, folder, network
            )
        )

    # Loop through each station, read in the file, gather the neccessary info, and append to the csv
    # Progress bar will be printed to console for each loop interation
    for i in progressbar(range(len(network_filenames))):
        filename = network_filenames[i]

        # Get full filepath
        filepath_s3 = "s3://{0}/{1}/{2}/{3}".format(bucket, folder, network, filename)

        # Read in the zarr
        ds = xr.open_zarr(filepath_s3)

        # Get lat and lon values from dataset
        # Confirm that there is only one value for lat and one for lon
        lon = np.unique(ds.lon.values)
        lat = np.unique(ds.lat.values)
        _warn_if_more_than_one_value_in_array(l=lon, name="longitude")
        _warn_if_more_than_one_value_in_array(l=lat, name="latitude")
        lon = lon[0]
        lat = lat[0]

        # Get elevation value
        # Try to grab from the dataset... if none found, set elevation to Nan
        # Confirm that only one value for elevation is found
        try:
            elevation = ds["elevation"].values
        except:
            elevation = np.nan
        if type(elevation) == np.ndarray:
            elevation = np.unique(elevation)
            _warn_if_more_than_one_value_in_array(l=elevation, name="elevation")
            elevation = elevation[0]

        # Append to dataframe as a new row
        df = pd.DataFrame(
            {
                "era-id": [filename.split(".zarr/")[0]],
                "elevation": [elevation],
                "network": [network],
                "longitude": [lon],
                "latitude": [lat],
            }
        )
        stations_df = pd.concat([stations_df, df], ignore_index=True)

    # Re-sort into alphabetical order, just in case!
    stations_df = stations_df.sort_values("era-id", ignore_index=True)

    # Rewrite/save locally
    stations_df.to_csv(csv_local_path, index=False)
    print("csv saved locally to: {0}".format(csv_local_path))

    # Rewrite/save to s3 bucket
    stations_df.to_csv(csv_s3_filepath, index=False)
    print("csv saved to s3 path: {0}".format(csv_s3_filepath))


if __name__ == "__main__":
    with warnings.catch_warnings():
        # Annoying RuntimeWarning is raised... ignore here so we can better evaluate the print output of the function, which includes a progress bar
        warnings.simplefilter("ignore")
        add_stations_to_csv()
