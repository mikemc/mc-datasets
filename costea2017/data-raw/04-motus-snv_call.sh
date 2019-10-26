#!/bin/bash
#SBATCH -o sb-motus-snv_call-%j.out
#SBATCH -c 4
set -euo pipefail

# Replace N with number of processors / cores you want to use
# If submitting to a cluster with sbatch, also adjust in "SBATCH -c N" above
nproc=4

# Set paths based on the DATA_PATH set in the .env file
export $(cat .env | xargs)
bam_path=$DATA_PATH/costea2017/motus/bam
out_path=$DATA_PATH/costea2017/motus/snv-calls

# Activate conda environment with motus. Must be edited based on system
# configuration
eval "$(conda shell.zsh hook)"
conda activate motus

# Call SNVs
motus snv_call \
    -d $bam_path \
    -o $out_path \
    -t $nproc
