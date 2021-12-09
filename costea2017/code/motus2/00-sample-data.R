# Setup -----------------------------------------------------------------------

library(tidyverse)
library(here)
# Also used:
# readxl
# dotenv

dotenv::load_dot_env(here(".env"))

# Location for downloading files and storing intermediate files. Create a
# symlink if want to store large data (the sequence reads) in a different
# location)
dl_path <- file.path(Sys.getenv("DATA_PATH"), "costea2017", "metadata")
if (!dir.exists(dl_path)) {
    dir.create(dl_path, recursive = TRUE)
}

# Download sample information -------------------------------------------------

## Download the sample metadata from the supplemental file provided by the
## authors as Excel spreadsheets
# Supplementary Data 2 : Members and composition of mock community
# https://media.nature.com/original/nature-assets/nbt/journal/v35/n11/extref/nbt.3960-S5.xlsx
# Supplementary Data 3 : Sample description
# https://media.nature.com/original/nature-assets/nbt/journal/v35/n11/extref/nbt.3960-S6.xlsx
urls <- c("https://media.nature.com/original/nature-assets/nbt/journal/v35/n11/extref/nbt.3960-S5.xlsx",
    "https://media.nature.com/original/nature-assets/nbt/journal/v35/n11/extref/nbt.3960-S6.xlsx")
names(urls) <- c("mock_composition.xlsx", "sample_description.xlsx")
fns <- file.path(dl_path, names(urls))
walk2(urls, fns, download.file)

## Sequence file metadata from the ENA
# https://www.ebi.ac.uk/ena/data/view/PRJEB14847
fn <- file.path(dl_path, "PRJEB14847.tsv")
download.file("https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJEB14847&result=read_run&download=txt", fn)

list.files(dl_path)

#  Mock community composition -------------------------------------------------
mock <- readxl::read_xlsx(file.path(dl_path, "mock_composition.xlsx"), 
    sheet = 1) %>%
    rename(Taxon = `Bacterial species`) %>%
    mutate(Taxon = str_extract(Taxon, "\\S+ \\S+"),
        Taxon = str_replace(Taxon, " ", "_")) %>%
    select(-`NCBI tax ID`)
head(mock)
costea2017_mock_composition <- mock
# TODO: save this somewhere

# Sample metadata from the manuscript supplementary files ---------------------

# The metadata for the three experimental phases of the Costea2017 experiment
# is given on a single sheet in separate rectagular ranges.
fn <- file.path(dl_path, "sample_description.xlsx")
ranges <- c("1" = "A3:B192", "2" = "D3:F77", "3" = "H3:J32")
sam.all <- map_dfr(ranges, ~readxl::read_xlsx(fn, range = .),
    .id = "phase") %>%
    rename_all(str_to_lower) %>%
    rename(si_sample = sample) %>%
    mutate(
        individual = case_when(
            phase == 3 ~ individual,
            TRUE ~ str_sub(si_sample, 1, 1)
            )
        ) %>%
    mutate_at(vars(phase, lab, protocol, individual), as.factor)
sam.all %>%
    group_by(phase) %>%
    count

# Link sample metadata to sequence data ---------------------------------------

# We have to use the tsv from the ENA to link the SI sample names to the ENA
# accessions
seqtb <- readr::read_tsv(file.path(dl_path, "PRJEB14847.tsv")) %>%
    select_if(~length(unique(.)) > 1) %>%
    select(-starts_with("sra"), -starts_with("fastq"), 
        -starts_with("submitted"))


# Note, there are a number of samples with multiple runs
seqtb %>%
    group_by(sample_accession) %>%
    count %>%
    filter(n>1)

a1 <- seqtb %>%
    group_by(sample_accession) %>%
    mutate(n = n()) %>%
    filter(n>1) %>%
    arrange(sample_accession, run_accession)

# GOAL: get the sample_accession matched to si_sample


# Sample names grouped by phase
si_sns <- sam.all %>%
    select(phase, si_sample) %>%
    group_by(phase) %>%
    nest %>%
    deframe %>%
    map(pull)

# si_sns$`1`
# si_sns$`2`



# The phase 3 si_sample can be found in: 
# library_name, experiment_alias, run_alias
# and has the pattern "BYQ_[A-Z]{4}"
# The phase 2 si_sample can be found at the end of the sample_alias col and has
# the pattern "(A|B)1_[0-9]{3}$"

tb <- seqtb %>%
    select(sample_accession, library_name, experiment_alias, run_alias,
        sample_alias) %>%
    mutate_all(str_replace_all, "[_ -]+", "_") %>%
    mutate(
        phase2_name = str_extract(sample_alias, "(A|B)1_[0-9]{3}$"),
        phase3_name = str_extract(library_name, "BYQ_[A-Z]{4}")
        ) %>%
    # Get the si sample name iff only one of the phase 2 or 3 names exists
    filter(!(is.na(phase2_name) & is.na(phase3_name))) %>%
    mutate(si_sample = case_when(
            !is.na(phase2_name) & is.na(phase3_name) ~ phase2_name,
            is.na(phase2_name) & !is.na(phase3_name) ~ phase3_name,
            TRUE ~ NA_character_
            )
        ) %>%
    select(sample_accession, si_sample) %>%
    # Remove duplicates
    distinct %>%
    # join with the sample data
    left_join(sam.all, by = "si_sample") %>%
    # Drop the remaining Phase 1 samples
    filter(phase %in% c("2", "3"))

# Are all phase 2 and 3 samples represented?
sam.all %>%
    filter(phase != 1) %>%
    pull(si_sample) %>%
    {all(. %in% tb$si_sample)}

# Check for multiple sample accessions per si_sample

tb %>%
    group_by(phase) %>%
    count

tb %>%
    group_by(phase, si_sample) %>%
    count %>%
    filter(n > 1)

tb %>%
    filter(is.na(phase))

# SAMEA4347651	B1_062	2	1	Q	B

a2 <- seqtb %>%
    # select(sample_accession, library_name, experiment_alias, run_alias,
    #     sample_alias) %>%
    mutate_all(str_replace_all, "[_ -]+", "_") %>%
    mutate(
        phase2_name = str_extract(sample_alias, "(A|B)1_[0-9]{3}$"),
        phase3_name = str_extract(library_name, "BYQ_[A-Z]{4}")
        ) %>%
    filter(phase2_name == "B1_062")
    
a4 <- seqtb %>%
    arrange(sample_alias) %>%
    group_by(sample_alias) %>%
    mutate(n = n()) %>%
    filter(n > 1)

# Oh. One of these is a phase 1 sample, B1_062_FNOSW .
# I have to reexamine my original method, and make sure I discard the phase 1
# samples 


# Q: how to distinguish phase 2 from phase 1?

# All of the library substrings for phase 1 have 5 characters
phase1_substrings <- sam.all %>%
    filter(phase == 1) %>%
    mutate(library_substring = str_extract(si_sample, "[A-Z]{4,6}$")) %>%
    pull(library_substring)
phase1_substrings %>%
    map_int(nchar) %>%
    summary
# The NA corresponds to a typo of 1 instead of I

# Get info from the library_name and sample_alias that can be used to filter
# out the "C" samples and the phase 1 samples, and get the si_sample name for
# phase 2 and phase 3 samples
# which phase
tb <- seqtb %>%
    select(sample_accession, library_name, sample_alias) %>%
    mutate_all(str_replace_all, "[_ -]+", "_") %>%
    mutate(
        library_substring = str_extract(library_name, 
            "(?<=AWF_)[A-Z]{4,}(?=_r)"),
        individual = str_extract(sample_alias, "(?<=_)(A|B|C)(?=(1|2)(_|$))"),
        phase2_name = str_extract(sample_alias, "(A|B)1_[0-9]{3}$"),
        phase3_name = str_extract(library_name, "BYQ_[A-Z]{4}")
        )
tb0 <- tb %>%
    # Remove Individual C samples
    filter((individual != "C") | is.na(individual)) %>%
    # Remove Phase 1 samples; necessary to avoid cross-matching with phase 2
    # names in sam.all
    filter(!(library_substring %in% phase1_substrings)) %>%
    # Pick the sample name via the phase2 or phase3 form
    mutate(si_sample = case_when(
            !is.na(phase2_name) & is.na(phase3_name) ~ phase2_name,
            is.na(phase2_name) & !is.na(phase3_name) ~ phase3_name,
            TRUE ~ NA_character_
            )
        ) %>%
    select(sample_accession, si_sample) %>%
    # Remove duplicates (check later)
    distinct %>%
    # join with the sample data
    left_join(sam.all, by = "si_sample")
tb0$si_sample %>% is.na %>% any
#> [1] FALSE

# Note, we have a 1-1 correspondence between sample accesions and si samples.
# but there can be multiple runs per sample.

# Ok, I think that tb0 is a trustworthy map from sample_accession to si_sample
# for phase 2 and 3 samples. (For phase 1 samples, more care is needed, and
# would have to fix the putative typos)

# Next: 
# - Call tb0 something sensible
# - compare tb0 w/ seqtb to see what's up with cases where multiple runs per
# si_sample // sample_accession

# One row per run, with relevant info
runtb <- seqtb %>%
    select(run_accession, sample_accession, instrument_model, nominal_length,
        library_layout, read_count, base_count, ) %>%
    left_join(tb0, by = "sample_accession") %>%
    filter(!is.na(phase)) %>%
    arrange(phase, lab, individual, protocol, sample_accession) %>%
    group_by(si_sample) %>%
    mutate(n = n()) %>%
    ungroup

# Duplicates - 
dups <- runtb %>%
    filter(n>1)

# Any difference in insert size or read count?
runtb %>%
    group_by(n, sample_accession) %>%
    summarize_at(vars(nominal_length, read_count), min) %>%
    summarize_at(vars(nominal_length, read_count), mean)
# In cases where there are two runs, the read counts tend to be lower; perhaps
# these samples were resequenced after the initial read count was found to be
# low?


# TODO: save runtb. Either use this for the metadata for the motus2
# output, or re-run motus on the duplicate samples with all the sequence data.

# revise the above; explain the reasoning

write_csv(runtb, here("costea2017", "data", "costea2017-run-metadata.csv"))
saveRDS(runtb, here("costea2017", "data", "costea2017-run-metadata.Rds"))

# OLD ----------



# Note, we expect 189 samples for phase 1 and 74 samples for phase 2,
# indicating extra samples or duplicate samples.

dups <- seqtb %>%
    filter(duplicated(SI_sample) | 
            duplicated(SI_sample, fromLast = TRUE)
        ) %>%
    arrange(Phase, SI_sample, run_accession)
# The library_name and sample_alias are duplicated; while the run_alias and
# run_accessions differ
dups %>%
    select(Phase, SI_sample, library_name, sample_alias) %>%
    group_by(Phase) %>%
    top_n(4, SI_sample)
dups %>%
    select(Phase, SI_sample, run_accession, run_alias) %>%
    group_by(Phase) %>%
    top_n(4, SI_sample)
# (Differences appear in the HAYK8ADXX part of the run alias)
# For Phase 1, there can be differences in library layout and length
dups %>%
    select(Phase, SI_sample, library_layout, nominal_length) %>%
    group_by(Phase) %>%
    top_n(4, SI_sample)

# Unclear why there are multiple libraries for these samples.

# Create final sample metadata table ------------------------------------------

# Let's get sample data for just Phases 2 and 3

# Get a minimal set of variables for joining with the sample metadata
seqtb0 <- seqtb %>%
    select(SI_sample, Phase, run_accession, 
        instrument_model, library_layout, nominal_length) %>%
    rename_at(vars(-SI_sample), str_to_sentence)
        # library_name, run_alias,
        # fastq_bytes, fastq_md5, fastq_ftp, fastq_aspera)

# Once we join with the sample data, we will have multiple rows for the
# SI_sample that are duplicated in the seqtb.
sam <- sam.all %>%
    filter(Phase %in% c(2, 3)) %>%
    left_join(seqtb0, by = c("Phase", "SI_sample")) %>%
    arrange(Phase, SI_sample, Run_accession) %>%
    group_by(SI_sample) %>%
    mutate(
        id = rank(Run_accession),
        n = n(),
        Sample = ifelse(n > 1,
            paste(SI_sample, id, sep = "_"), 
            SI_sample),
        ) %>%
    select(Sample, Phase, everything())

sam %>%
    # filter(duplicated(SI_sample) | duplicated(SI_sample, fromLast = TRUE)) %>%
    filter(n > 1) %>%
    select(Phase, Sample, SI_sample, Run_accession, id, n)

sam <- sam %>%
    mutate(
        Individual = case_when(
            Phase == 2 ~ str_sub(SI_sample, 1, 1),
            Phase == 3 ~ Individual),
        Sample = case_when(
            # Phase == 2 ~ paste0("Ph2", "L", Lab, Protocol, Individual),
            Phase == 2 ~ Sample,
            Phase == 3 ~ paste0(Protocol, Individual))
        )

write_csv(sam, here("costea2017", "data", "costea2017-phases23-sample-data.csv"))




####################### OLD below here; should save ####################







# Link sample metadata to sequence data ---------------------------------------

# We have to use the tsv from the ENA to link the SI sample names to the ENA
# accessions
seqtb <- readr::read_tsv(file.path(dl_path, "PRJEB14847.tsv")) %>%
    select_if(~length(unique(.)) > 1) %>%
    select(-starts_with("sra"), -starts_with("fastq"), 
        -starts_with("submitted"))

# There are a few strings for each run that will help us in matching with the
# metadata
seqtb <- seqtb %>%
    mutate(
        String1 = str_extract(library_name, "(AWF)|(BYQ)"),
        String2 = str_extract(library_name, "(?<=(AWF|BYQ)_)[A-Z]{4,6}"),
        String3 = str_extract(sample_alias, 
            "([ABC][1I][_ -]{1,2}[:digit:]{3}([_ -][:digit:])?)|C[12]"),
    )
# seqtb$String1
# seqtb$String2
# seqtb$String3

# Note, there is considerable heterogeneity in the formatting of String3,
# including an apparent typo of AI that should be A1. The first letter [ABC]
# indicates (I think) the specimen. That is followed by a 1, an I (typo), or a
# 2, with the 2 only appearing for C2. 

# We are not interested in the C samples since these are not included in the
# sample metadata so let's drop these now. 
seqtb <- seqtb %>%
    filter(is.na(String3) | (str_sub(String3, 1, 1) != "C"))


# For the purposes of matching to the samples, let's fix up String3: fix the
# typo and make the separator formatting (of underscore, space, and dash)
# uniform 
seqtb <- seqtb %>%
    mutate(
        String3 = str_replace(String3, "AI", "A1"),
        String3 = str_replace_all(String3, "[_ -]+", "_"),
    )

# Challenge is to use the Strings to identify the phase

# Phase 3 samples are easy to identify and get the sample names of
seqtb <- seqtb %>%
    mutate(
        Phase = ifelse(String1 == "BYQ", "3", NA),
        SI_sample = ifelse(Phase == "3", paste("BYQ", String2, sep = "_"), NA),
    )
seqtb %>%
    filter(Phase == 3) %>%
    .$SI_sample

# But for Phase 1 and 2 we can't easily tell, and will need to match against
# the sample metadata. 

## To do this, first create the names we would see if these samples were from
## the particular phase. 
seqtb <- seqtb %>%
    mutate(
        Name_if_1 = paste(String3, String2, sep = "_"),
        Name_if_2 = String3,
    )
## Then, look for matches in the known sample names.
samples.phase1 <- sam.all %>% filter(Phase == 1) %>% .$SI_sample
samples.phase2 <- sam.all %>% filter(Phase == 2) %>% .$SI_sample
# Check what's missing in Phase 1
setdiff(samples.phase1, seqtb$Name_if_1)
# Both appear to be typos. We can find these two missing samples if we suppose:
# - Sample A1_081_A1OSW should be A1_081_AIOSW (typo of 1 <-> I)
# - Sample A1_054_HLOSW corresponds to String2 = HLOSW & String3 = A1_045; a
#   typo of (054 <- 045).

# Let's adjust Name_if_1 to match the sample metadata
seqtb <- seqtb %>%
    mutate(Name_if_1 = case_when(
        Name_if_1 == "A1_081_AIOSW" ~ "A1_081_A1OSW",
        Name_if_1 == "A1_045_HLOSW" ~ "A1_054_HLOSW",
        TRUE ~ Name_if_1
        )
    )
setdiff(samples.phase1, seqtb$Name_if_1)

# Check what's missing in Phase 2
setdiff(samples.phase2, seqtb$Name_if_2)
# all good!

seqtb <- seqtb %>%
    mutate(
        Phase = case_when(
            Name_if_1 %in% samples.phase1 ~ "1",
            Name_if_2 %in% samples.phase2 ~ "2",
            TRUE ~ Phase),
        SI_sample = case_when(
            Phase == 1 ~ Name_if_1,
            Phase == 2 ~ Name_if_2,
            Phase == 3 ~ SI_sample)
        ) %>%
    arrange(Phase)
# Check that the extracted names look good
seqtb %>%
    group_by(Phase) %>%
    top_n(5, SI_sample) %>%
    select(Phase, SI_sample)
# Check the counts by Phase and make sure no missing sample names
seqtb %>%
    group_by(Phase, is.na(SI_sample)) %>%
    count

# Note, we expect 189 samples for phase 1 and 74 samples for phase 2,
# indicating extra samples or duplicate samples.

dups <- seqtb %>%
    filter(duplicated(SI_sample) | 
            duplicated(SI_sample, fromLast = TRUE)
        ) %>%
    arrange(Phase, SI_sample, run_accession)
# The library_name and sample_alias are duplicated; while the run_alias and
# run_accessions differ
dups %>%
    select(Phase, SI_sample, library_name, sample_alias) %>%
    group_by(Phase) %>%
    top_n(4, SI_sample)
dups %>%
    select(Phase, SI_sample, run_accession, run_alias) %>%
    group_by(Phase) %>%
    top_n(4, SI_sample)
# (Differences appear in the HAYK8ADXX part of the run alias)
# For Phase 1, there can be differences in library layout and length
dups %>%
    select(Phase, SI_sample, library_layout, nominal_length) %>%
    group_by(Phase) %>%
    top_n(4, SI_sample)

# Unclear why there are multiple libraries for these samples.

# Create final sample metadata table ------------------------------------------

# Let's get sample data for just Phases 2 and 3

# Get a minimal set of variables for joining with the sample metadata
seqtb0 <- seqtb %>%
    select(SI_sample, Phase, run_accession, 
        instrument_model, library_layout, nominal_length) %>%
    rename_at(vars(-SI_sample), str_to_sentence)
        # library_name, run_alias,
        # fastq_bytes, fastq_md5, fastq_ftp, fastq_aspera)

# Once we join with the sample data, we will have multiple rows for the
# SI_sample that are duplicated in the seqtb.
sam <- sam.all %>%
    filter(Phase %in% c(2, 3)) %>%
    left_join(seqtb0, by = c("Phase", "SI_sample")) %>%
    arrange(Phase, SI_sample, Run_accession) %>%
    group_by(SI_sample) %>%
    mutate(
        id = rank(Run_accession),
        n = n(),
        Sample = ifelse(n > 1,
            paste(SI_sample, id, sep = "_"), 
            SI_sample),
        ) %>%
    select(Sample, Phase, everything())

sam %>%
    # filter(duplicated(SI_sample) | duplicated(SI_sample, fromLast = TRUE)) %>%
    filter(n > 1) %>%
    select(Phase, Sample, SI_sample, Run_accession, id, n)

sam <- sam %>%
    mutate(
        Individual = case_when(
            Phase == 2 ~ str_sub(SI_sample, 1, 1),
            Phase == 3 ~ Individual),
        Sample = case_when(
            # Phase == 2 ~ paste0("Ph2", "L", Lab, Protocol, Individual),
            Phase == 2 ~ Sample,
            Phase == 3 ~ paste0(Protocol, Individual))
        )

write_csv(sam, here("costea2017", "data", "costea2017-phases23-sample-data.csv"))

