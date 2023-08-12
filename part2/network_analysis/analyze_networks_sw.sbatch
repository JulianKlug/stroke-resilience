#!/bin/bash
#SBATCH --partition=public-longrun-cpu
#SBATCH --job-name=resilience_analyze_networks_sw    # Job name
#SBATCH --ntasks=1
#SBATCH --output=/home/users/k/klug/logs/resilience/compute_sw_%j.log   # Standard output and error log
#SBATCH --licenses=matlab@matlablm.unige.ch
#SBATCH --time=10-00:00:00               # Time limit hrs:min:sec


module load MATLAB/2022a

BASE_MFILE_NAME='analyze_networks_sw'
MFILE_DIR='/home/users/k/klug/stroke-resilience/part2/network_analysis/'

cd ${MFILE_DIR}

unset DISPLAY
echo "Running ${BASE_MFILE_NAME}.m on $(hostname)"
srun matlab -nodesktop -nosplash -nodisplay -r ${BASE_MFILE_NAME}
echo "Finished running ${BASE_MFILE_NAME}.m"

