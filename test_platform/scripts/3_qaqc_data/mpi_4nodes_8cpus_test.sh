#!/bin/bash -l

# Job Information:
# Name: hist-obs MPI Job
# Description: This script runs an MPI-optimized job using Conda environment "hist-obs" on a cluster with 4 nodes and 8 CPUs.
# Partition: test (for 4 nodes, 8 CPUs)
# Time Limit: 5 hours
# Python Script: ALLNETWORKS_qaqc.py
# Arguments: network="CAHYDRO"

#SBATCH --job-name=hist-obs          # Name of the job
#SBATCH --ntasks=8                   # Total number of MPI processes (1 per CPU)
#SBATCH --nodes=4                    # Number of nodes (1 node with 2 CPUs per node)
#SBATCH --ntasks-per-node=2          # Number of MPI processes per node (2 processes per node)
#SBATCH --time=5:00:00               # Maximum runtime (adjust as necessary)
#SBATCH --partition=test             # Queue/partition (adjust as necessary)
#SBATCH --output=%x_%j_output.txt    # Standard output file with job name and job ID
#SBATCH --error=%x_%j_error.txt      # Standard error file with job name and job ID

# AWS secret info 
export AWS_ACCESS_KEY_ID="put-your-key-id-here"
export AWS_SECRET_ACCESS_KEY="put-your-key-here"
export AWS_DEFAULT_REGION="us-west-2"

# Load OpenMPI
module load openmpi

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

# Start time
start_time=$(date +%s)

# Run the Python script directly using Python from the activated environment
srun --mpi=pmi2 python3 ${PYSCRIPT} --network="CAHYDRO" 

# End time (right after the job finishes)
end_time=$(date +%s)

# Calculate elapsed time
elapsed_time=$((end_time - start_time))
echo "Job completed in $elapsed_time seconds."

# Deactivate Conda environment after job completion
conda deactivate