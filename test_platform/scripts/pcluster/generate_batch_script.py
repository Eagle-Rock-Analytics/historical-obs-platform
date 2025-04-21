"""
generate_batch_script.py

Reads a SLURM batch script template and fills in the network name and number
of stations to set up an embarrassingly parallel array job. Splits large arrays
into chunks of 1000 tasks and generates multiple SLURM batch scripts accordingly.

Example usage:
---------------
python generate_batch_script.py --network=LOXWFO

This will read the template `run_qaqc_template.sh`, replace placeholders 
({NETWORK}, {ARRAY_RANGE}), and write one or more batch scripts like:
- `run_qaqc_LOXWFO.sh`       (if ≤1000 stations)
- `run_qaqc_LOXWFO_1.sh`     (if >1000, multiple scripts)
"""

import argparse
from pathlib import Path
import math


def main():
    parser = argparse.ArgumentParser(
        description="Generate SLURM batch script(s) for a network, splitting by 1000 tasks if needed."
    )
    parser.add_argument(
        "-n", "--network", required=True, help="Network name (e.g., LOXWFO)"
    )
    args = parser.parse_args()
    network = args.network

    # Paths
    input_file = Path(f"stations_input/{network}-input.dat")
    template_file = Path("run_qaqc_template.sh")

    if not input_file.exists():
        raise FileNotFoundError(f"Missing input file: {input_file}")
    if not template_file.exists():
        raise FileNotFoundError(f"Missing template file: {template_file}")

    # Count stations
    with input_file.open("r") as f:
        nrows = sum(1 for _ in f)

    # Read template
    with template_file.open("r") as f:
        template = f.read()

    # Determine chunking
    chunk_size = 1000
    if nrows <= chunk_size:
        array_range = f"1-{nrows}"
        script = template.replace("{NETWORK}", network).replace(
            "{ARRAY_RANGE}", array_range
        )
        output_file = Path(f"run_qaqc_{network}.sh")
        with output_file.open("w") as f:
            f.write(script)
        print(f"Batch script written to {output_file}")
    else:
        num_chunks = math.ceil(nrows / chunk_size)
        for i in range(num_chunks):
            start = i * chunk_size + 1
            end = min((i + 1) * chunk_size, nrows)
            array_range = f"{start}-{end}"
            script = template.replace("{NETWORK}", network).replace(
                "{ARRAY_RANGE}", array_range
            )
            output_file = Path(f"run_qaqc_{network}_{i+1}.sh")
            with output_file.open("w") as f:
                f.write(script)
            print(f"Batch script written to {output_file}")


if __name__ == "__main__":
    main()
