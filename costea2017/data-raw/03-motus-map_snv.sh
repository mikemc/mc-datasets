#!/bin/bash
#SBATCH -o sb-motus-map_snv-%j.out
#SBATCH -c 4
set -euo pipefail

# Replace N with number of processors / cores you want to use. If submitting
# to a cluster with sbatch, also adjust in "SBATCH -c N" above
nproc=4

# Command-line args:
# Path with the reads
reads_path=$1
# Path for bam output 
out_path=$2
# The ERA run accession, used to create the fastq.gz file names
accession=$3

# Activate conda environment with motus. Must be edited based on system
# configuration
eval "$(conda shell.zsh hook)"
conda activate motus

# Align reads
motus map_snv \
    -f $reads_path/${accession}_1.fastq.gz \
    -r $reads_path/${accession}_2.fastq.gz \
    -t $nproc \
    > $out_path/${accession}.bam
