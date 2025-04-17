"""
generate_station_list.py

This script generates a list of station IDs for a specified network and saves them to a .dat file.

Overview:
---------
1. Filters stations from a CSV file based on the specified network.
2. Extracts the 'era-id' values for the filtered stations.
3. Writes the station IDs to a .dat file, one per line, without quotes.

Inputs:
-------
- Network name (string) to filter the station data.
- CSV file containing station data (must include 'network' and 'era-id' columns).

Outputs:
--------
- A .dat file containing the list of station IDs for the specified network.

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
    Generates a list of station IDs for a specified network and saves it to a `.dat` file.

    This script filters stations from a cleaned CSV dataset based on the provided network and
    extracts their corresponding 'era-id' values. The results are then written to a `.dat` file
    with each station ID on a new line.

    Parameters:
    -----------
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
    network_df = stations_df[stations_df["network"] == network]

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

    # Define the file path for the output .dat file (saved to the current directory)
    filename = f"{stations_input_dir}/{network}-input.dat"

    # Write the station IDs to the file, one per line, without quotes
    with open(filename, "w") as f:
        for i, era_id in enumerate(era_ids):
            f.write(str(era_id))
            if i < len(era_ids) - 1:  # Add newline only if it's not the last item
                f.write("\n")

    print(f"Station list successfully written to {filename}")

    print("\nUpdate the following lines in the batch script:")
    print(f"#SBATCH --array=1-{num_stations}")
    print(
        f'STATION=$(awk "NR==$SLURM_ARRAY_TASK_ID" stations_input/{network}-input.dat)\n'
    )

    return None


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
