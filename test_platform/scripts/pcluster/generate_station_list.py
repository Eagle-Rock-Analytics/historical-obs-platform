"""
generate_station_list.py

This script generates a list of station IDs for a specified network and saves them to one or more .dat files.

Overview:
---------
1. Filters stations from a CSV file based on the specified network.
2. Extracts the 'era-id' values for the filtered stations.
3. Writes the station IDs to one or more .dat files, one per line, without quotes.
   If the number of stations exceeds 1000, multiple files will be created.

Inputs:
-------
- Network name (string) to filter the station data.
- CSV file containing station data (must include 'network' and 'era-id' columns).

Outputs:
--------
- One or more .dat files containing the list of station IDs for the specified network.

Example usage:
--------------
python generate_station_list.py -n "LOXWFO"
python generate_station_list.py --network=LOXWFO
"""

import pandas as pd
import argparse
from pathlib import Path


def generate_station_list(network: str):
    """
    Generates a list of station IDs for a specified network and saves to `.dat` files.

    If there are more than 1000 stations, the list is split into multiple files with
    names like `{network}_1-input.dat`, `{network}_2-input.dat`, etc.

    Parameters
    ----------
    network : str
        The network to process. Corresponds to the network name in the dataset.
    """
    # Read the CSV file containing station data
    csv_filepath = "s3://wecc-historical-wx/2_clean_wx/temp_clean_all_station_list.csv"
    stations_df = pd.read_csv(csv_filepath)

    # Define and create the directory (if it doesn't already exist)
    stations_input_dir = Path("stations_input")
    stations_input_dir.mkdir(parents=True, exist_ok=True)

    # Filter the dataframe to only include rows corresponding to the specified network
    # And, only cleaned stations
    network_df = stations_df[
        (stations_df["network"] == network) & (stations_df["cleaned"] == "Y")
    ]

    # Check if nothing is returned. Raise ValueError and print useful message.
    if len(network_df) == 0:
        unique_networks = ", ".join(stations_df["network"].unique())  # Unique networks
        raise ValueError(
            f"No stations found for network: {network}. Available networks: {unique_networks}"
        )

    # Get the 'era-id' column as a list (array of station IDs)
    era_ids = sorted(
        network_df["era-id"].values
    )  # Sort into alphabetical/numerical order
    num_stations = len(era_ids)
    print(f"{num_stations} stations found.")

    # If fewer than or equal to 1000 stations, save to a single file
    chunk_size = 1000
    file_summary = []

    if num_stations <= chunk_size:
        filename = stations_input_dir / f"{network}-input.dat"

        # Write the station IDs to the file, one per line, without quotes
        with open(filename, "w") as f:
            f.write("\n".join(map(str, era_ids)))

        print(f"Station list successfully written to {filename}")
        file_summary.append((filename.name, 1, num_stations))

    # If more than 1000 stations, split into chunks
    else:
        for i in range(0, num_stations, chunk_size):
            chunk = era_ids[i:i + chunk_size]
            index = i // chunk_size + 1
            start_idx = i + 1
            end_idx = min(i + chunk_size, num_stations)
            filename = stations_input_dir / f"{network}_{index}-input.dat"

            # Write the station IDs to the file
            with open(filename, "w") as f:
                f.write("\n".join(map(str, chunk)))

            print(f"Chunk {index} written to {filename}")
            file_summary.append((filename.name, start_idx, end_idx))

    # Print a summary of all the files written
    print("\nSummary of station input files:")
    for fname, start, end in file_summary:
        print(f"  - {fname}: stations {start}–{end} ({end - start + 1} total)")
    print(f"\n{len(file_summary)} file(s) written.")


def main():
    """
    The main function to parse arguments and call `generate_station_list`.

    Command-line arguments:
    ------------------------
    -n, --network    : The network to process (e.g., 'LOXWFO')
    """
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="Generate a station list for a given network"
    )

    # Define arguments (network is required, input CSV is hardcoded)
    parser.add_argument(
        "-n",
        "--network",
        required=True,
        help="Network name to filter stations (e.g., 'LOXWFO').",
    )

    # Parse arguments
    args = parser.parse_args()

    # Call the function to generate the station list
    print(f"Generating station input list for network: {args.network}")
    generate_station_list(args.network)
    print("Script complete.")


if __name__ == "__main__":
    main()