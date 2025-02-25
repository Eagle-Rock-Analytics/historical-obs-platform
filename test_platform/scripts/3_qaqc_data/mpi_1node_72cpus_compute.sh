#!/bin/bash -l

# Job Information:
# Name: hist-obs MPI Job
# Description: This script runs an MPI-optimized job using Conda environment "hist-obs" on a single node with 72 CPUs.
# Partition: compute (for 1 node, 72 CPUs)
# Time Limit: 5 hours
# Python Script: ALLNETWORKS_qaqc.py
# Arguments: network="CAHYDRO"
# 
# NOTE: This script must be run from the directory containing ALLNETWORKS_qaqc.py.

#SBATCH --job-name=hist-obs            # Name of the job
#SBATCH --ntasks=72                    # Total number of MPI processes (1 per CPU)
#SBATCH --nodes=1                      # Number of nodes (1 node with 72 CPUs)
#SBATCH --ntasks-per-node=72           # Number of MPI processes per node (72 processes per node)
#SBATCH --time=5:00:00                 # Maximum runtime (adjust as necessary)
#SBATCH --partition=compute            # Queue/partition (adjust as necessary)
#SBATCH --output=%x_%j_output.txt      # Standard output file with job name and job ID
#SBATCH --error=%x_%j_error.txt        # Standard error file with job name and job ID

# AWS secret info 
export AWS_ACCESS_KEY_ID="put-your-key-id-here"
export AWS_SECRET_ACCESS_KEY="put-your-key-here"
export AWS_DEFAULT_REGION="us-west-2"

# Load OpenMPI
module load openmpi

# Print SBATCH job settings for debugging
echo "====================================="
echo "Job Name: $SLURM_JOB_NAME"
echo "Job ID: $SLURM_JOB_ID"
echo "Partition: $SLURM_JOB_PARTITION"
echo "Number of Nodes: $SLURM_JOB_NUM_NODES"
echo "Tasks Per Node: $SLURM_NTASKS_PER_NODE"
echo "Total Tasks: $SLURM_NTASKS"
echo "CPUs Per Task: $SLURM_CPUS_PER_TASK"
echo "Job Start Time: $(date)"
echo "====================================="

# Define the path to your Python script
PYSCRIPT="ALLNETWORKS_qaqc.py"

# Check if the Python script exists
if [ ! -f "$PYSCRIPT" ]; then
  echo "Error: Python script not found!"
  exit 1
fi

# Load Conda initialization 
source /shared/miniconda3/etc/profile.d/conda.sh

# Activate the Conda environment
conda activate /shared/miniconda3/envs/hist-obs

# Start time tracking
start_time=$(date +%s)

# Run the Python script directly using Python from the activated environment
srun --mpi=pmi2 python3 ${PYSCRIPT} --network="CAHYDRO" 

# End time tracking
end_time=$(date +%s)

# Calculate and print elapsed time
elapsed_time=$((end_time - start_time))
echo "====================================="
echo "Job completed in $elapsed_time seconds."
echo "Job End Time: $(date)"
echo "====================================="

# Deactivate Conda environment after job completion
conda deactivate