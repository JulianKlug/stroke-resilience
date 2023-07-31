#!/bin/bash
#SBATCH --job-name=resilience_analyze_networks_sw    # Job name
#SBATCH --ntasks=1
#SBATCH --mem=100G                     # Job memory request
#SBATCH --cpus-per-task=16
#SBATCH --output=/home/users/k/klug/logs/resilience/compute_sw_%j.log   # Standard output and error log
#SBATCH --licenses=matlab@matlablm.unige.ch


module load MATLAB/2022a

BASE_MFILE_NAME='/home/users/k/klug/stroke-resilience/part2/network_analysis/analyze_networks_sw.m'

unset DISPLAY

echo "Running ${BASE_MFILE_NAME}.m on $(hostname)"

srun matlab -nodesktop -nosplash -nodisplay -r ${BASE_MFILE_NAME}
