#!/bin/bash -l

# Job Information:
# Name: hist-obs MPI Job
# Description: This script runs an MPI-optimized job using Conda environment "hist-obs" on a cluster with 4 nodes and 8 CPUs.
# Partition: test (for 4 nodes, 8 CPUs)
# Time Limit: 5 hours
# Python Script: /shared/nicole/historical-obs-platform/test_platform/scripts/3_qaqc_data/ALLNETWORKS_qaqc.py
# Arguments: network="CAHYDRO"

#SBATCH --job-name=hist-obs          # Name of the job
#SBATCH --ntasks=8                   # Total number of MPI processes (1 per CPU)
#SBATCH --nodes=4                    # Number of nodes (1 node with 2 CPUs per node)
#SBATCH --ntasks-per-node=2          # Number of MPI processes per node (2 processes per node)
#SBATCH --time=5:00:00               # Maximum runtime (adjust as necessary)
#SBATCH --partition=test             # Queue/partition (adjust as necessary)
#SBATCH --output=mpi_job_output.txt  # Standard output file
#SBATCH --error=mpi_job_error.txt    # Standard error file

# Load required modules (OpenMPI and Conda)
module load openmpi

# Activate conda env 
conda activate hist-obs

# Define the path to your Python script
PYSCRIPT=/shared/nicole/historical-obs-platform/test_platform/scripts/3_qaqc_data/ALLNETWORKS_qaqc.py

# Run the Python script using Conda environment
srun --mpi=pmix_v3 conda run -n hist-obs python3 ${PYSCRIPT} network="CAHYDRO"

# Deactivate Conda environment after job completion
conda deactivate