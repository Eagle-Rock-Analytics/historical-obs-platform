#!/bin/bash -l

#!/bin/bash -l

# Job Information:
# Name: hist-obs MPI Job
# Description: This script runs an MPI-optimized job using Conda environment "hist-obs" on a cluster with 4 nodes and 8 CPUs.
# Partition: test (for 4 nodes, 8 CPUs)
# Time Limit: 5 hours
# Python Script: /shared/nicole/historical-obs-platform/test_platform/scripts/3_qaqc_data/ALLNETWORKS_qaqc.py
# Arguments: network="CAHYDRO"

# Set output and error file paths
OUTPUT_FILE="${HOME}/%x_%j_output.txt"
ERROR_FILE="${HOME}/%x_%j_error.txt"

#SBATCH --job-name=hist-obs          # Name of the job
#SBATCH --ntasks=8                   # Total number of MPI processes (1 per CPU)
#SBATCH --nodes=4                    # Number of nodes (1 node with 2 CPUs per node)
#SBATCH --ntasks-per-node=2          # Number of MPI processes per node (2 processes per node)
#SBATCH --time=5:00:00               # Maximum runtime (adjust as necessary)
#SBATCH --partition=test             # Queue/partition (adjust as necessary)
#SBATCH --output=$OUTPUT_FILE        # Standard output file with job name and job ID
#SBATCH --error=$ERROR_FILE          # Standard error file with job name and job ID

# Load required modules (OpenMPI and Conda)
module load openmpi

# Initialize Conda
source /shared/miniconda3/etc/profile.d/conda.sh   # Correct path to Conda initialization

# Activate Conda environment (full path to environment)
conda activate /shared/miniconda3/envs/hist-obs   # Ensure no extra "/envs" is added

# Define the path to your Python script
PYSCRIPT="${HOME}/historical-obs-platform/test_platform/scripts/3_qaqc_data/ALLNETWORKS_qaqc.py"

# Run the Python script using Conda environment with MPI plugin version 5 
srun --mpi=pmix_v5 conda run -n hist-obs python3 ${PYSCRIPT} network="CAHYDRO"

# Deactivate Conda environment after job completion
conda deactivate