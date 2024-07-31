#!/bin/bash -l

#SBATCH -N 16
#SBATCH -t 06:00:00
#SBATCH -o stdout/training_qaqc_3.out
#SBATCH -e err/training_qaqc_3.err
#SBATCH -p compute-16vcores
#SBATCH -J train_qaqc
#SBATCH --ntasks=256

echo ""
echo "STARTING PREAMBLE"
date
echo ""

# Start measuring time
start_time=$(date +%s.%N)

# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/shared/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/shared/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/shared/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/shared/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

sleep 2;
conda activate /shared/anaconda3/envs/hist-obs
sleep 2;
export PATH="/shared/anaconda3/envs/hist-obs/bin":$PATH

cd /shared/historical-obs-platform/test_platform/scripts/3_qaqc_data/

echo ""
echo "STARTING CALCULATION"
date
echo ""

mpirun -n 288 python3 Training_QAQC.py

echo ""
echo "CALCULATION ENDED"
date
echo ""

# End measuring time
end_time=$(date +%s.%N)

# Calculate elapsed time
elapsed_time=$(echo "$end_time - $start_time" | bc)
elapsed_time=$(echo "$elapsed_time / 1" | bc)
# Calculate hours
hours=`printf "%02d" "$(( ${elapsed_time} / 3600 ))"`
remainder=`printf "%02d" "$(( ${elapsed_time} % 3600 ))"`
# Calculate minutes
minutes=`printf "%02d" "$(( remainder / 60 ))"`
# Calculate seconds
seconds=`printf "%02d" "$(( remainder % 60 ))"`

echo "Elapsed time: ${hours}:${minutes}:${seconds}"

sleep 5;
exit
