"""
generate_batch_script.py

Reads a SLURM batch script template (run_{process}_template.sh) and fills in
the network name and number of stations to set up an embarrassingly parallel
array job. The --process argument selects which template and output filename
to use:

  - "qaqc"  -> template: run_qaqc_template.sh,  output: run_qaqc_{network}.sh
  - "merge"  -> template: run_merge_template.sh, output: run_merge_{network}.sh

Placeholders replaced in the template:
  {NETWORK} -> network name
  {NROWS}   -> number of stations in stations_input/{network}-input.dat

Example usage:
---------------
python generate_batch_script.py --network=LOXWFO --process=qaqc
python generate_batch_script.py --network=LOXWFO --process=merge
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
    parser.add_argument(
        "-p",
        "--process",
        required=True,
        choices=["merge", "qaqc"],
        help="Process type: 'merge' or 'qaqc'",
    )
    args = parser.parse_args()
    network = args.network
    process = args.process

    # Paths
    input_file = Path(f"stations_input/{network}-input.dat")
    template_file = Path(f"run_{process}_template.sh")
    output_file = Path(f"run_{process}_{network}.sh")

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
