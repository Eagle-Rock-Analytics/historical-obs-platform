#!/bin/bash -l

################################################################################
# SLURM Batch Script Template: run_qaqc_missing.sh
#
# Description:
#   Launches the QA/QC pipeline for historical weather station data.
#   Each SLURM array task processes a single station using the script:
#       ../3_qaqc_data/QAQC_run_for_single_station.py
#
#   This is an embarrassingly parallel workload:
#   - Each task is independent and processes one station
#   - Tasks require no communication or coordination
#   - Ideal for SLURM array jobs and horizontal scaling
#
# Inputs:
#   - Station list: stations_input/missing-input.dat
#   - Python script: ../3_qaqc_data/QAQC_run_for_single_station.py
#   - Conda environment: hist-obs
#
# SLURM Configuration:
#   - Partition: compute-72cpus
#   - CPUs per task: 1
#   - One array task per station (set using line count of input file)
#
# Working Directory:
#   This script should be submitted from the `pcluster/` directory.
#
# Output:
#   - Output and error logs are saved per-task, including the station name
#
# Notes:
#   - AWS credentials must be exported before submission (or hardcoded below)
#   - The array range will be set automatically based on station count
################################################################################

# Job Information:
#SBATCH --job-name=hist-obs
#SBATCH --array=1-22
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --partition=compute
#SBATCH --output=%x_%A_%a_output.txt
#SBATCH --error=%x_%A_%a_error.txt

# Get the station name for this array task
STATION=$(awk "NR==$SLURM_ARRAY_TASK_ID" stations_input/missing-input.dat)

# AWS credentials
# Don't need to hard code them in if they are already saved as environment variables
# export AWS_ACCESS_KEY_ID="put-your-key-id-here"
# export AWS_SECRET_ACCESS_KEY="put-your-key-here"
# export AWS_DEFAULT_REGION="us-west-2"

# Rename SLURM-generated output and error files to include station name
ORIG_OUT="${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_output.txt"
ORIG_ERR="${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_error.txt"
NEW_OUT="${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${STATION}_output.txt"
NEW_ERR="${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${STATION}_error.txt"

# Wait for SLURM to generate the default files, then rename
sleep 2
[ -f "$ORIG_OUT" ] && mv "$ORIG_OUT" "$NEW_OUT"
[ -f "$ORIG_ERR" ] && mv "$ORIG_ERR" "$NEW_ERR"

# Use new output filename as log file
log_file="$NEW_OUT"

# Print SBATCH job settings for debugging to the log file
{
  echo "====================================="
  echo "Station: $STATION"
  echo "Job Name: $SLURM_JOB_NAME"
  echo "Array Job ID: $SLURM_ARRAY_JOB_ID"
  echo "Task ID: $SLURM_ARRAY_TASK_ID"
  echo "Partition: $SLURM_JOB_PARTITION"
  echo "Number of Nodes: $SLURM_JOB_NUM_NODES"
  echo "Tasks Per Node: $SLURM_NTASKS_PER_NODE"
  echo "Total Tasks: $SLURM_NTASKS"
  echo "CPUs Per Task: $SLURM_CPUS_PER_TASK"
  echo "Job Start Time: $(date)"
  echo "====================================="
} >> "$log_file"

# Change to the directory containing the script
cd ../3_qaqc_data/ || { echo "Directory change failed"; exit 1; }

# Define the path to your Python script
PYSCRIPT="QAQC_run_for_single_station.py"

# Start time tracking
start_time=$(date +%s)

# Load Conda initialization
source /shared/nicole/.mamba/etc/profile.d/conda.sh

# Run the Python script
source /opt/parallelcluster/shared/miniforge3/b/etc/profile.d/conda.sh
conda activate /shared/nicole/.mamba/envs/hist-obs
python3 ${PYSCRIPT} --station="$STATION"

# End time tracking
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

# Change to the directory containing the logfile
cd ../pcluster/ || { echo "Directory change failed"; exit 1; }

# Write end-of-job info
{
  echo "====================================="
  echo "Job completed in $elapsed_time seconds."
  echo "Job End Time: $(date)"
  echo "====================================="
} >> "$log_file"
