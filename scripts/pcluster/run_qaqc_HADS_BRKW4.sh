#!/bin/bash -l

################################################################################
# SLURM Batch Script: run_qaqc_HADS_BRKW4.sh
#
# Description:
#   Reruns the QA/QC pipeline for stations that failed due to
#   deaccumulation plot crash: HADS_BRKW4, CWOP_E6545
################################################################################

# Job Information:
#SBATCH --job-name=hist-obs
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --partition=compute
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

# AWS credentials
# Don't need to hard code them in if they are already saved as environment variables
# export AWS_ACCESS_KEY_ID="put-your-key-id-here"
# export AWS_SECRET_ACCESS_KEY="put-your-key-here"
# export AWS_DEFAULT_REGION="us-west-2"

log_file="${SLURM_JOB_NAME}_${SLURM_JOB_ID}_output.txt"

# Print SBATCH job settings for debugging to the log file
{
  echo "====================================="
  echo "Job Name: $SLURM_JOB_NAME"
  echo "Job ID: $SLURM_JOB_ID"
  echo "Partition: $SLURM_JOB_PARTITION"
  echo "Job Start Time: $(date)"
  echo "====================================="
} >> "$log_file"

# Change to the directory containing the script
cd ../3_qaqc_data/ || { echo "Directory change failed"; exit 1; }

# Define the path to your Python script
PYSCRIPT="QAQC_run_for_single_station.py"

# Load Conda initialization
source /shared/nicole/.mamba/etc/profile.d/conda.sh
source /opt/parallelcluster/shared/miniforge3/b/etc/profile.d/conda.sh
conda activate /shared/nicole/.mamba/envs/hist-obs

# Run stations in serial
for STATION in HADS_BRKW4 CWOP_E6545; do
    echo "====================================="
    echo "Running QAQC for: $STATION"
    echo "Start Time: $(date)"
    echo "====================================="

    start_time=$(date +%s)
    python3 ${PYSCRIPT} --station="$STATION"
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))

    echo "Completed $STATION in $elapsed_time seconds."
    echo ""
done

# Change to the directory containing the logfile
cd ../pcluster/ || { echo "Directory change failed"; exit 1; }

# Write end-of-job info
{
  echo "====================================="
  echo "All stations completed."
  echo "Job End Time: $(date)"
  echo "====================================="
} >> "$log_file"
