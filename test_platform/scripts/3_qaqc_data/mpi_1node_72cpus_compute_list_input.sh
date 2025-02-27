#!/bin/bash -l

# Job Information:
# Name: hist-obs MPI Job
# Description: This script runs an MPI-optimized job using Conda environment "hist-obs" on a single node with 72 CPUs.
# Partition: compute (for 1 node, 72 CPUs)
# Time Limit: 5 hours
# Python Script: ALLNETWORKS_qaqc.py
# Arguments: Each network from networks-input.dat is processed in a separate job array task.
#
# NOTE: This script must be run from the directory containing ALLNETWORKS_qaqc.py.

#SBATCH --job-name=hist-obs            # Name of the job
#SBATCH --array=1-$(wc -l < networks-input.dat)  # Create a job array for each network
#SBATCH --ntasks=1                     # Each array job runs a single task
#SBATCH --cpus-per-task=72             # Number of CPUs per task
#SBATCH --time=5:00:00                 # Maximum runtime (adjust as necessary)
#SBATCH --partition=compute            # Queue/partition (adjust as necessary)
#SBATCH --output=stdout/%x_%A_%a.out   # Standard output file per array task
#SBATCH --error=err/%x_%A_%a.err       # Standard error file per array task

# AWS secret info 
export AWS_ACCESS_KEY_ID="put-your-key-id-here"  # AWS access key
export AWS_SECRET_ACCESS_KEY="put-your-key-here"  # AWS secret key
export AWS_DEFAULT_REGION="us-west-2"  # AWS region

# Load OpenMPI for parallel processing
module load openmpi

# Define the path to your Python script
PYSCRIPT="ALLNETWORKS_qaqc.py"

# Get the network name for this array task
NETWORK=$(sed -n "${SLURM_ARRAY_TASK_ID}p" networks-input.dat)

# Create a log file based on job name, job ID, and task ID
log_file="${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_output.txt"

# Print SBATCH job settings for debugging to the log file
{
  echo "====================================="
  echo "Job Name: $SLURM_JOB_NAME"
  echo "Job ID: $SLURM_JOB_ID"
  echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
  echo "Processing Network: $NETWORK"
  echo "Partition: $SLURM_JOB_PARTITION"
  echo "CPUs Per Task: $SLURM_CPUS_PER_TASK"
  echo "Job Start Time: $(date)"
  echo "====================================="
} >> "$log_file"

# Check if the Python script exists
if [ ! -f "$PYSCRIPT" ]; then
  echo "Error: Python script not found!" >> "$log_file"
  exit 1
fi

# Load Conda initialization 
source /shared/miniconda3/etc/profile.d/conda.sh

# Activate the Conda environment
conda activate /shared/miniconda3/envs/hist-obs

# Start time tracking
start_time=$(date +%s)

# Run the Python script with the selected network
srun --mpi=pmi2 python3 ${PYSCRIPT} --network=${NETWORK} --sample=1

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

# Deactivate Conda environment after job completion
conda deactivate