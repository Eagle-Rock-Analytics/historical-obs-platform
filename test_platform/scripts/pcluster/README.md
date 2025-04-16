# Running the QAQC Script on the Cluster (with Slurm)
Here's a useful guide for how to run the QAQC code in an AWS pcluster environment  🖥️✨

## 📖 Read This First!

- **Changing Stations?**  
  See the [Changing the Stations to Run](#changing-the-stations-to-run) section below.
  
- **PLEASE always remember to turn off the compute fleet after you've finished your jobs!**   
  We get charged by AWS when the fleet is on. 
---

## 🚀 Step-by-Step Instructions

1. **Turn on the compute fleet**  

2. **Log in to the cluster:**  
```
ssh hist-obs-cluster
```
3. **Navigate to the scripts directory:**  
```
cd historical-obs-platform/test_platform/scripts/pcluster
```

4. **Update the batch script:**  
Make sure your AWS credentials included! But, ⚠️ Do **NOT** push this modified batch script to GitHub since it contains your private info! You'll need to modify the following rows in the batch script with your unique info: 
```
export AWS_ACCESS_KEY_ID="put-your-key-id-here"
export AWS_SECRET_ACCESS_KEY="put-your-key-here"
export AWS_DEFAULT_REGION="us-west-2"
```   
5. **Submit the batch script:**  
```
sbatch run.sh
```
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
