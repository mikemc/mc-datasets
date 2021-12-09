# Submit `kneaddata` jobs to Slurm --------------------------------------------

library(tidyverse)
library(here)
library(fs)

# Load sample metadata with run accessions
# sam <- here("costea2017", "output", "sample-data", 
#   "costea2017-phases23-sample-data.csv"
# ) %>%
#   read_csv %>%
#   janitor::clean_names()

#### Set up input and output files

raw_dir <- here("costea2017/data/reads/raw")
out_dir <- here("costea2017/data/reads/kneaddata")
db_dir <- here("data/kneaddata/Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1")
dir_create(out_dir)

x <- tibble(path = dir_ls(raw_dir, glob = '*.fastq.gz')) %>%
  mutate(
    file = path_file(path),
    sample_id_dir = str_extract(file, '^.+(?=\\.fastq\\.gz$)'),
  ) %>%
  separate(sample_id_dir, into = c('sample_id', 'direction'), sep = '_') %>%
  select(-file) %>%
  pivot_wider(names_from = direction, names_prefix = 'r', values_from = path) %>%
  glimpse

#### Generate sbatch commands

# Example kneaddata command:
# $ kneaddata --input1 seq1.fastq --input2 seq2.fastq -db $DATABASE --output kneaddata_output
# Note, fastqc does not detect adapters in the test file.

cmds <- x %>%
  mutate(
    knead_cmd = str_glue(
      'kneaddata --input1 {r1} --input2 {r2} -db {db_dir} --output {out_dir}'
    ),
    job_name = str_glue('{sample_id}-kneaddata'),
    sb_cmd = str_glue(
      'sbatch --exclusive --mem=0 --job-name={job_name} --output=slurm-%j-{job_name}.out --wrap="{knead_cmd}"'
    )
  ) %>%
  pull(sb_cmd)

#### Submit jobs

walk(commands, system)
