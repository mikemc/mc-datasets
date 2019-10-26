# Submit `motus map_snv` jobs to Slurm ----------------------------------------

library(readr)
library(here)

# library(dotenv)
# load_dot_env(here(".env"))
# data_path <- Sys.getenv("DATA_PATH")

# Load sample metadata with run accessions
sam <- here("costea2017", "data", "costea2017-phases23-sample-data.csv") %>%
    read_csv

# Generate commands for calling mOTUs2 map_snv
sam <- sam %>%
    mutate(command = paste("sbatch 03-motus-map_snv.sh", Run_accession))

walk(sam$command, system)
