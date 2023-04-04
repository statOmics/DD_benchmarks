## Filter lupus data samples: replicates and underrepresented categories
## =====================================================================


# Setup -------------------------------------------------------------------

## Directory setup
here_root <- "benchmarks/lupus"
here::i_am(file.path(here_root, "scripts", "01-filter-samples.R"))

data_dir <- here::here(here_root, "data")
lupus_file <- file.path(data_dir, "lupus-SCE-cleaned.rds")
stopifnot(file.exists(lupus_file))

out_file <- file.path(data_dir, "lupus-SCE-filtered.rds")


## Libraries
suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(HDF5Array)
    library(scuttle)
    library(BiocParallel)
})

## Load data
lupus <- readRDS(lupus_file)


# Sample filtering --------------------------------------------------------

## Load sample metadata
sample_metadata <- metadata(lupus)$sample_metadata

## Remove underrepresented categories
sample_metadata <- subset(sample_metadata,
    Sex == "Female" & Ethnicity %in% c("Asian", "European") & Status != "Flare"
)
sample_metadata <- droplevels(sample_metadata)
keep_samples <- as.character(sample_metadata$sample)

## Subset and update SCE object
lupus_filtered <- lupus[, lupus$sample %in% keep_samples]
colData(lupus_filtered) <- droplevels(colData(lupus_filtered))

## Add categorical age variable
sample_age <- sample_metadata$Age
sample_age <- as.numeric(levels(sample_age))[sample_age]

age_breaks <- c(20, 40, 60, 80)
age_cat <- cut(sample_age, age_breaks, right = FALSE)

sample_metadata$Age <- sample_age
sample_metadata$age_cat <- age_cat

lupus_filtered$age_cat <- sample_metadata[lupus_filtered$sample, "age_cat"]

metadata(lupus_filtered)$sample_metadata <- sample_metadata


# Split patients ----------------------------------------------------------

## Split patients per Ethnicity
patients <- as.character(sample_metadata$ind_cov)
patient_splits <- split(patients, sample_metadata$Ethnicity)
metadata(lupus_filtered)$patient_splits <- patient_splits


# Export output data ------------------------------------------------------

saveRDS(lupus_filtered, file = out_file)
