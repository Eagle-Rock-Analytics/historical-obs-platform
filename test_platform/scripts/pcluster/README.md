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

### To change stations
 1. Create a file named `{NETWORK}-input.dat` containing the stations you wish to run. You can use the provided Python script `generate_station_list.py` to generate this file—refer to the script header for detailed instructions. Each station should be listed on a new line, without commas or quotes, and **must** be stations in the same network. Before creating the file, check if one already exists for the network you're interested in, as some have been pre-generated 😊.
      Example: 
      ```   
      LOXWFO_CRXC1   
      LOXWFO_ABC12
      ```   
  2. Update the batch script’s array size to match the number of rows:   
       ```bash
       #SBATCH --array=1-n
       ```
       Replace `n` with the total number of stations (rows). For example, if you have two stations, use `1-2`.
  3. Update the filename in the batch script to match `{NETWORK}-input.dat`. For example, for LOXWFO:
      ```bash
      STATIONS_INPUT="LOXWFO-input.dat"
      ```

### Important 
 - The station string must match the actual station ID in AWS.
 - No quotation marks around network names (e.g., `LOXWFO_CRXC1`, not `"LOXWFO_CRXC1"`).

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
