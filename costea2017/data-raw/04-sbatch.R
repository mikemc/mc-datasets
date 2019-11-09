# Submit `motus snv_call` jobs to Slurm ----------------------------------------

# library(readr)
# library(dplyr)
# library(purrr)
library(here)
library(dotenv)

load_dot_env(here(".env"))
data_path <- Sys.getenv("DATA_PATH")

# # Load sample metadata with run accessions
# sam <- here("costea2017", "data", "costea2017-phases23-sample-data.csv") %>%
#     read_csv

script_path <- here("costea2017", "data-raw", "04-motus-snv_call.sh")
bam_path <- file.path(data_path, "costea2017/motus/bam")
out_path <- file.path(data_path, "costea2017/motus/2019-10-28-snv_call")

# TODO: make sure out path doesn't already exist; consider giving a date-stamp on it as well
# TODO: consider if this script is necessary, since we're only submitting 1 job
# NOTE: this step only takes a few minutes

# Generate command for calling job script
command <- paste("sbatch", script_path, bam_path, out_path)

system(command)
