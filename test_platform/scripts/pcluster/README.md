# Running the QA/QC Script on the Cluster (with SLURM)

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
   Open `run_qaqc_array.sh` and update these lines:
   ```bash
   export AWS_ACCESS_KEY_ID="your-key-id"
   export AWS_SECRET_ACCESS_KEY="your-secret-key"
   export AWS_DEFAULT_REGION="us-west-2"
   But, ⚠️ Do **NOT** push this modified batch script to GitHub since it contains your private info!

5. **Submit the batch script to SLURM:**  
   Once your batch script is updated with the correct AWS credentials and station information, submit it using the following command:
   ```bash
   sbatch run_qaqc_array.sh

6. **Monitor your run** 
Use "squeue", or check the output/error files (see below).

7. **Turn off the compute fleet**  
(Don’t forget! 💡)

---
## 🔄 Changing the Stations to Run

- The batch script (`run.sh`) reads station info from `stations-input.dat`.
- **To change stations:**
1. Add the desired station(s) to `stations-input.dat`.  
  - One station per line, no commas or quotes.  
  - Example:  
    ```
    LOXWFO_CRXC1
    SFOXYZ_ABC12
    ```
2. Update the batch script’s array size to match the number of rows:  
  - Edit line 14:  
    ```
    #SBATCH --array=1-n
    ```
    Replace `n` with the total number of stations (rows).
  - For two networks, use `1-2`.

- **Important:**  
- The station string must match the actual station ID in AWS.
- No quotation marks around network names (e.g., `LOXWFO_CRXC1`, not `"LOXWFO_CRXC1"`).

---

## 📝 Reference Info

- **Check job duration:**  
- Look in the output file, usually in the same directory as the batch script.
- Output file format:  
 ```
 ${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${STATION}_output.txt
 ```
- Error file format:  
 ```
 ${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${STATION}_error.txt
 ```

- **Monitor jobs:**  
- Use `squeue` to see job status.

---

## 🌟 Tips & Reminders

- Double-check your AWS info before submitting!
- Don’t push sensitive info to GitHub.
- Keep your `stations-input.dat` and batch script in sync.
- If you get stuck, check the output/error files for clues! 🕵️‍♂️

---
