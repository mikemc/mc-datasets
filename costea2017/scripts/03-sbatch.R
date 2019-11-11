# Submit `motus map_snv` jobs to Slurm ----------------------------------------

library(readr)
library(dplyr)
library(purrr)
library(here)
library(dotenv)

load_dot_env(here(".env"))
data_path <- Sys.getenv("DATA_PATH")

# Load sample metadata with run accessions
sam <- here("costea2017", "data", "costea2017-phases23-sample-data.csv") %>%
    read_csv

script_path <- here("costea2017", "data-raw", "03-motus-map_snv.sh")
reads_path <- file.path(data_path, "costea2017/reads")
out_path <- file.path(data_path, "costea2017/motus/bam")

# Generate commands for calling mOTUs2 map_snv
commands <- paste("sbatch", script_path, reads_path, out_path,
    sam$Run_accession)

walk(commands, system)
