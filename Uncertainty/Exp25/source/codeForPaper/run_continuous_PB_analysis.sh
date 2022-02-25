#!/bin/sh
#SBATCH -p short
#SBATCH -t 3:00:00
#SBATCH --mem=128G
#SBATCH -c 6
#SBATCH --mail-user=melaniebasnak@g.harvard.edu

module load matlab/2018b
matlab -nodesktop -r "addpath('/n/scratch3/users/m/mb491/codes/'); continuous_PB_data_analysis('/n/scratch3/users/m/mb491/data/${1}',${2})"
