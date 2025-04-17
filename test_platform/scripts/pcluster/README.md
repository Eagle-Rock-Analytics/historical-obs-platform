# Running the QA/QC script on an AWS pcluster with SLURM

Welcome! This guide walks you through running the QA/QC pipeline for historical weather station data using the AWS ParallelCluster environment. 🛰️🖥️

---

## 📖 Read This First!

- **Want to process different stations?**  
  See the [🔄 Changing the Stations to Run](#changing-the-stations-to-run) section below.

- **PLEASE remember to turn off the compute fleet after your jobs finish!**  
  AWS charges while the fleet is running ⏳💸

---

## 🚀 Step-by-Step Instructions

1. **Start the compute fleet**  
   Use the AWS console or CLI to turn it on.

2. **SSH into the cluster:**  
   Use your unique cluster alias, which may or may not be the same as `hist-obs-cluster` depending on what you've set it to.  
   ```bash
   ssh hist-obs-cluster

3. **Navigate to the SLURM batch script directory:**  
   ```bash
   cd historical-obs-platform/test_platform/scripts/pcluster

4. **Edit the batch script with your AWS credentials:**  
   But, ⚠️ Do **NOT** push this modified batch script to GitHub since it contains your private info!   
   Open `run_qaqc_array.sh` and update these lines:
   ```bash
   export AWS_ACCESS_KEY_ID="your-key-id"
   export AWS_SECRET_ACCESS_KEY="your-secret-key"
   export AWS_DEFAULT_REGION="us-west-2"

6. **Submit the batch script to SLURM:**  
   Once your batch script is updated with the correct AWS credentials and station information, submit it using the following command:
   ```bash
   sbatch run_qaqc_array.sh

7. **Monitor your run** 
 - Use squeue to check the job queue
 - Or tail the log files in the logs/ folder).

8. **Turn off the compute fleet**  
💡 Please don't forget! This helps us avoid surprise charges. 

---
## 🔄 Changing the Stations to Run
Each SLURM array task processes one station in a fully independent and serial job — a classic [embarrassingly parallel](https://en.wikipedia.org/wiki/Embarrassingly_parallel) workload. The batch script (`run_qaqc_array.sh`) reads station info from `{NETWORK}-input.dat`, where {NETWORK} is the name of the weather station network (i.e. for network LOXWFO, the file is named `LOXWFO-input.dat`).

1. **Create a station list file**  
   Create a file named `{NETWORK}-input.dat` containing the stations you want to run. Each station should be on a new line—no commas, quotes, or extra formatting—and **all stations must belong to the same network**.

   👉 You can use the provided `generate_station_list.py` script to create this file. Check the script's header for usage instructions.

   Before creating a new file, check the `stations_input/` folder—some lists have already been pre-generated 😊.

   **Example file contents:**
   ```
   LOXWFO_CRXC1
   LOXWFO_ABC12
   ```

2. **Update the batch script**  
   The batch script needs two updates:
   - The array size (to match the number of stations)  
   - The path to the correct `{NETWORK}-input.dat` file

   You can do this in one of two ways:

   **Option 1 (Recommended): Use the generator script**
   ```bash
   python generate_batch_script -n="LOXWFO"
   ```
   This will automatically read from `LOXWFO-input.dat` and update the batch script for you.

   **Option 2: Manual edits**
   - Set the array size:
     ```bash
     #SBATCH --array=1-n
     ```
     Replace `n` with the number of stations (e.g., `1-2` if you have two stations).

   - Update the station file path:
     ```bash
     STATION=$(awk "NR==$SLURM_ARRAY_TASK_ID" stations_input/LOXWFO-input.dat)
     ```

### Important 
 - The station string must match the actual station ID in AWS.
 - No quotation marks around station ids in the `{NETWORK}-input.dat` (e.g., `LOXWFO_CRXC1`, not `"LOXWFO_CRXC1"`).
 - Each station id must be a new row in the `{NETWORK}-input.dat` file. 

---

## 📝 Reference Info

### Check job duration 
 - Look in the output file, usually in the same directory as the batch script.
 - Output file format:  
    ```
    ${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${STATION}_output.txt
    ```
 - Error file format:  
    ```
    ${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${STATION}_error.txt
    ```

### Monitor jobs 
 - Use `squeue` to see job status.

---

## 🌟 Tips & Reminders
✅ Double-check AWS credentials before submitting  
🔒 Don’t commit modified batch scripts with secrets to GitHub  
🔄 Keep your station list and SLURM array size in sync  
🧪 Use squeue and log files to debug if things go wrong  

---
