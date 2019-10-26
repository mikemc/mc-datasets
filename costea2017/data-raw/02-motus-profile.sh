#!/bin/bash
#SBATCH -o sb-motus-profile-%j.out
#SBATCH -c 4
set -euo pipefail

# Replace N with number of processors / cores you want to use. If submitting
# to a cluster with sbatch, also adjust in "SBATCH -c N" above
nproc=4

# Set paths based on the DATA_PATH set in the .env file
export $(cat .env | xargs)
reads_path=$DATA_PATH/costea2017/reads
out_path=$DATA_PATH/costea2017/motus

# The ERA run accession, used to create the fastq.gz file names; read in as
# command line arg
accession=$1

# Activate conda environment with motus. Must be edited based on system
# configuration
eval "$(conda shell.zsh hook)"
conda activate motus

# Run motus with "profile" option
motus profile \
    -f $reads_path/${accession}_1.fastq.gz \
    -r $reads_path/${accession}_2.fastq.gz \
    -o $out_path/${accession}-counts.motus \
    -M $out_path/${accession}.mgc \
    -n $accession \
    -c \
    -g 1 \
    -y insert.raw_counts \
    -t $nproc

# Afterwards, merge all with
#> FLIST=$(echo ERR*-counts.motus | tr ' ' ,)
#> motus merge -i $FLIST -o costea2017-counts.motus
