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

# Define the path to your Python script
PYSCRIPT="${HOME}/historical-obs-platform/test_platform/scripts/3_qaqc_data/ALLNETWORKS_qaqc.py"

# Check the current working directory
echo "Current working directory: $(pwd)"

# Check if the Python script exists
PYSCRIPT=/home/nicole/historical-obs-platform/test_platform/scripts/3_qaqc_data/ALLNETWORKS_qaqc.py
echo "Checking if Python script exists: $PYSCRIPT"
if [ ! -f "$PYSCRIPT" ]; then
  echo "Error: Python script not found!"
  exit 1
fi

# Load Conda initialization (if needed)
source /shared/miniconda3/etc/profile.d/conda.sh

# Activate the Conda environment (only once)
conda activate /shared/miniconda3/envs/hist-obs

# Run the Python script directly using Python from the activated environment
srun --mpi=pmix_v5 python3 ${PYSCRIPT} network="CAHYDRO"

# Deactivate Conda environment after job completion
conda deactivate