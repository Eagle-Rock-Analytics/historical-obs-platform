#!/bin/bash

# ==============================================================================
# Job Information:
# Name: mpi_local_4cpus_qaqc.sh
# Description: Runs the QA/QC process using OpenMPI on a local machine with 4 CPUs.
# Python Script: ALLNETWORKS_qaqc.py
# Arguments: --network="CAHYDRO"  # Example argument for network
# Note: This script utilizes 4 CPU cores (1 node with 4 CPUs) for parallel execution.
# To run this script, activate the conda environment, then run: 
#     bash mpi_1node_4cpus_local.sh
# 
# To modify the number of processors:
# - Change the value of "-np 4" in the line:
#     mpirun -np 4 python3 ${PYSCRIPT} --network="CAHYDRO"
# ==============================================================================

# Set up the log file with a timestamp
log_file="hist_obs_job_$(date +%Y%m%d_%H%M%S)_output.txt"

# Print script settings for debugging to the log file
{
  echo "====================================="
  echo "Job Start Time: $(date)"
  echo "Number of Processors: 4"
  echo "====================================="
} >> "$log_file"

# Define the path to your Python script
PYSCRIPT="ALLNETWORKS_qaqc.py"

# Check if the Python script exists
if [ ! -f "$PYSCRIPT" ]; then
  echo "Error: Python script not found!" >> "$log_file"
  exit 1
fi

# Start time tracking
start_time=$(date +%s)

# Run the Python script with OpenMPI
# Using 4 processes
# Modify the number of processors depending on your machine: see comment at top of file
echo "Running Python script with OpenMPI..." >> "$log_file"
mpirun -np 4 python3 ${PYSCRIPT} --network="CAHYDRO"

# End time tracking
end_time=$(date +%s)

# Calculate and print elapsed time to the log file
elapsed_time=$((end_time - start_time))
{
  echo "====================================="
  echo "Job completed in $elapsed_time seconds."
  echo "Job End Time: $(date)"
  echo "====================================="
} >> "$log_file"