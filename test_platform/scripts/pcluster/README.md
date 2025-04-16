# Historical Weather Station QA/QC Pipeline  
*SLURM Batch Processing Guide*

Welcome! This guide explains how to run the QA/QC pipeline for historical weather station data using AWS ParallelCluster and SLURM. It covers setup, job submission, and monitoring for efficient, parallelized processing.

---

## Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Step-by-Step Instructions](#step-by-step-instructions)
- [Changing the Stations to Run](#changing-the-stations-to-run)
- [Job Output and Monitoring](#job-output-and-monitoring)
- [Tips & Reminders](#tips--reminders)

---

## Overview

This pipeline processes historical weather station data by running a Python QA/QC script for each station in parallel, using a SLURM array job. Each SLURM task processes one station, making the workflow highly scalable and efficient.

---

## Prerequisites

- Access to an AWS ParallelCluster environment with SLURM.
- Your AWS credentials (do **not** share or commit these).
- A conda environment named `hist-obs` with required dependencies.
- The following directory structure (relative to your SSH login directory):

historical-obs-platform/
└── test_platform/
├── scripts/
│ └── pcluster/
│ └── run_qaqc_array.sh
├── 3_qaqc_data/
│ └── QAQC_run_for_single_station.py
└── stations_input/
└── stations-input.dat

text

---

## Step-by-Step Instructions

1. **Start the Compute Fleet**  
 Use the AWS Console or CLI to start the compute fleet.

2. **SSH into the Cluster**  
 Replace `hist-obs-cluster` with your cluster's alias:
ssh hist-obs-cluster

text

3. **Navigate to the Batch Script Directory**  
cd historical-obs-platform/test_platform/scripts/pcluster

text

4. **Edit the Batch Script with Your AWS Credentials**  
Open `run_qaqc_array.sh` and update these lines with your credentials:
export AWS_ACCESS_KEY_ID="your-key-id"
export AWS_SECRET_ACCESS_KEY="your-secret-key"
export AWS_DEFAULT_REGION="us-west-2"

text
**⚠️ Never commit or push this file with credentials to GitHub.**

5. **Submit the Batch Script to SLURM**  
Ensure your AWS credentials and station list are set, then run:
sbatch run_qaqc_array.sh

text

6. **Monitor Your Run**  
- Use `squeue` to check job status.
- Output and error logs are in the `logs/` directory.

7. **Turn Off the Compute Fleet**  
When jobs finish, stop the fleet to avoid unnecessary AWS charges.

---

## Changing the Stations to Run

- The batch script reads station IDs from `stations_input/stations-input.dat`.
- **To change stations:**
- Edit `stations-input.dat`, placing one station ID per line (no quotes or commas).
 ```
 LOXWFO_CRXC1
 SFOXYZ_ABC12
 ```
- Count the number of stations (lines) and update the array size in `run_qaqc_array.sh`:
 ```
 #SBATCH --array=1-N
 ```
 Replace `N` with the total number of stations.

- **Important:**  
- Station IDs must match those in AWS.
- Do **not** use quotation marks or extra spaces.

---

## Job Output and Monitoring

- **Logs:**  
Output and error logs are saved in the `logs/` directory with filenames:
{JOB_NAME}{ARRAY_JOB_ID}{TASK_ID}{STATION}output.txt
{JOB_NAME}{ARRAY_JOB_ID}{TASK_ID}_{STATION}_error.txt

text
- **Check Job Duration:**  
Open the output log file to view job start/end times and elapsed time.

- **Monitor Jobs:**  
Use:
squeue

text
to view job status.

---

## Tips & Reminders

- Double-check AWS credentials before submitting.
- Never push sensitive info (like AWS keys) to GitHub.
- Keep `stations-input.dat` and `--array` setting in sync.
- If errors occur, check the corresponding output/error files in `logs/`.
- Always turn off the compute fleet when done.

---

*For more on README formatting and best practices, see [GitHub Docs] and [README Best Practices].*

---

**Current date:** Wednesday, April 16, 2025

---

Let me know if you need a section on troubleshooting or further customization!