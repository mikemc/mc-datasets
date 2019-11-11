#!/bin/bash
#SBATCH -o sb-motus-snv_call-%j.out
#SBATCH -c 8
set -euo pipefail

# Replace N with number of processors / cores you want to use
# If submitting to a cluster with sbatch, also adjust in "SBATCH -c N" above
nproc=8

# Command-line args:
# Path with bam files produced by mapping step
bam_path=$1
# Path for snv-calling output 
out_path=$2

# Activate conda environment with motus. Must be edited based on system
# configuration
eval "$(conda shell.zsh hook)"
conda activate motus

# Call SNVs
# motus snv_call \
#     -d $bam_path \
#     -o $out_path \
#     -t $nproc
motus snv_call \
    -d $bam_path \
    -o $out_path \
    -fb 20 \
    -fd 1 \
    -fp 0.2 \
    -fc 1 \
    -t $nproc
