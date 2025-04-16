#!/bin/bash -l

# Job Information:
# Name: hist-obs MPI Job
# Description: This script runs an MPI-optimized job using Conda environment "hist-obs" on a single node with 72 CPUs.
# Partition: compute (for 1 node, 72 CPUs)
# Time Limit: 72 hours
# Python Script: QAQC_run_for_single_station.py
# Arguments: station=[row from stations-input.dat]

# NOTE: This script must be submitted from the directory containing stations-input.dat

# Job Information:
#SBATCH --job-name=hist-obs            # Name of the job
#SBATCH --array=1-4                    # Number of array jobs 
#SBATCH --ntasks=1                     # One task per array job
#SBATCH --nodes=1                      # All tasks run on 1 node
#SBATCH --ntasks-per-node=1            # 1 task per CPU 
#SBATCH --time=2:00:00                 # Time limit for the job
#SBATCH --partition=compute-72cpus     # Queue/partition to use
#SBATCH --output=logs/%x_%A_%a_output.txt   # Output file for each array job (task ID %a will be replaced by the array index)
#SBATCH --error=logs/%x_%A_%a_error.txt     # Error file for each array job

# Get the station name for this array task from the stations_input folder
STATIONS_INPUT="stations-input.dat"
STATION=$(awk "NR==$SLURM_ARRAY_TASK_ID" "stations_input/$STATIONS_INPUT")

# AWS credentials
export AWS_ACCESS_KEY_ID="put-your-key-id-here"
export AWS_SECRET_ACCESS_KEY="put-your-key-here"
export AWS_DEFAULT_REGION="us-west-2"

# Create the logs folder if it doesn't exist
mkdir -p logs

# Rename SLURM-generated output and error files to include station name
ORIG_OUT="logs/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_output.txt"
ORIG_ERR="logs/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_error.txt"
NEW_OUT="logs/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${STATION}_output.txt"
NEW_ERR="logs/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${STATION}_error.txt"

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

# Name of python script 
PYSCRIPT="QAQC_run_for_single_station.py"

# Start time tracking
start_time=$(date +%s)

# Run the Python script using conda
conda run -p /shared/miniconda3/envs/hist-obs python3 ${PYSCRIPT} --station="$STATION"

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