# Download reads from the ENA for Phase 2 and Phase 3 -------------------------

library(tidyverse)
library(here)
library(dotenv)

## Setup

# Load .env file
load_dot_env(here(".env"))
# Directory for downloading project data
data_path <- Sys.getenv("DATA_PATH")
# Paths to the aspera connect `ascp` program and the private key file to be
# used with the -i option
Sys.getenv("ASPERA_ASCP")
Sys.getenv("ASPERA_KEY")

# TODO: consider putting all the intermediate files into DATA_PATH as well

reads_path <- file.path(data_path, "costea2017", "reads")
if (!dir.exists(reads_path)) {
    dir.create(reads_path)
}

# Set data frame with aspera links for accessions for phase 2 and 3 data
sam <- here("costea2017", "data", "costea2017-phases23-sample-data.csv") %>%
    read_csv %>%
    rename(run_accession = Run_accession)
seqtb <- file.path(data_path, "costea2017", "metadata", "PRJEB14847.tsv") %>%
    read_tsv %>%
    left_join(sam, by = "run_accession") %>%
    filter(Phase %in% c(2, 3))

# TODO: edit the above once I make the sample data colnames lowercase

## Download with ascp (aspera connect command line tool)
# If don't have aspera, can download with wget (see below).
# The aspera urls are in the format "url/for/read1;url/for/read2" in `seqtb`
# and so we first split all out into a single list
# aspera_urls <- seqtb$fastq_aspera %>% str_split(";", simplify=TRUE) %>% c
aspera_urls <- seqtb %>%
    pull(fastq_aspera) %>% 
    str_split(";", simplify=TRUE) %>%
    c
commands <- paste(
    Sys.getenv("ASPERA_ASCP"),
    "-QT -l 300m -P33001 -i", 
    Sys.getenv("ASPERA_KEY"),
    paste0("era-fasp@", aspera_urls),
    reads_path
    )
walk(commands, system)

# ## Alternately, download with wget:
# ftp_urls <- seqtb$fastq_ftp %>% str_split(";", simplify=TRUE) %>% c
# dir.create(file.path(data_path, "reads"), recursive = TRUE)
# commands <- paste("wget", "-P", file.path(data_path, "reads"), 
#     paste0("ftp://", ftp_urls)
# walk(commands, system)

## Check dowloaded files against md5sums
downloads <- aspera_urls %>% 
    str_extract("ERR[0-9]*_[1-2]\\.fastq\\.gz") %>%
    file.path(reads_path, .)
md5sums_expected <- seqtb$fastq_md5 %>% str_split(";", simplify=TRUE) %>% c
md5sums_actual <- tools::md5sum(downloads)
# Fraction of files that were successfully downloaded and gave an md5
mean(!is.na(md5sums_actual))
# Check that the downloaded files match the expected md5
all(md5sums_expected == md5sums_actual, na.rm = TRUE)
