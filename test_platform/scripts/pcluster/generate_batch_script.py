"""
generate_batch_script.py

Reads a SLURM batch script template and fills in the network name and number
of stations to set up an embarrassingly parallel array job. Updates the header
comment to reflect the specific network.

Usage:
    python generate_batch_script.py --network NETWORK_NAME

This will read the template `run_qaqc_template.sh`, replace placeholders 
({network}, {nrows}), and write a batch script `run_qaqc_{network}.sh`.
"""

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Generate SLURM batch script for a network."
    )
    parser.add_argument(
        "-n", "--network", required=True, help="Network name (e.g., LOXWFO)"
    )
    args = parser.parse_args()
    network = args.network

    # Paths
    input_file = Path(f"stations_input/{network}-input.dat")
    template_file = Path("run_qaqc_template.sh")
    output_file = Path(f"run_qaqc_{network}.sh")

    if not input_file.exists():
        raise FileNotFoundError(f"Missing input file: {input_file}")
    if not template_file.exists():
        raise FileNotFoundError(f"Missing template file: {template_file}")

    # Count stations
    with input_file.open("r") as f:
        nrows = sum(1 for _ in f)

    # Read and modify template
    with template_file.open("r") as f:
        template = f.read()

    # Replace {NETWORK} with network name
    # Replace {NROWS} with string number of rows
    script = template.replace("{NETWORK}", network).replace("{NROWS}", str(nrows))

    # Write final script
    with output_file.open("w") as f:
        f.write(script)

    print(f"Batch script written to {output_file}")


if __name__ == "__main__":
    main()
