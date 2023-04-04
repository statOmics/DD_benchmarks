## Prepare input data for mock analyses and celltype subsets
## =========================================================

# Command line arguments --------------------------------------------------

library(argparse)
parser <- ArgumentParser()
parser$add_argument("-v", "--verbose",
    action = "store_true",
    help = "Print verbose output. Optional flag."
)
parser$add_argument("-i", "--infile",
    type = "character",
    default = "data/lupus-SCE-filtered.rds",
    help = "Input data file name. Default: \'%(default)s\'."
)
parser$add_argument("-o", "--outfile",
    type = "character",
    default = "data/lupus-SCE-mock-h5/se.rds",
    help = "Output file name for the SummarizedExperiment. Default: \'%(default)s\'."
)

args <- parser$parse_args()
verbose <- args$verbose

if (verbose) {
    message("Using input file: ", args$infile)
    message("Using output file: ", args$outfile)
}

## Directory setup
here_root <- "benchmarks/lupus"
here::i_am(file.path(here_root, "scripts", "02-prepare-mock.R"))

in_file <- here::here(here_root, args$infile)
h5_out_dir <- here::here(here_root, dirname(args$outfile))

if (!file.exists(in_file)) stop("Could not find input file ", in_file)
stopifnot(dir.exists(dirname(h5_out_dir)))

## Libraries
suppressPackageStartupMessages({
    library(HDF5Array)
    library(scuttle)
    library(dplyr)
})

if (verbose) message("Loading data...")
lupus <- readRDS(in_file)

seed <- 20220505


# Select samples ----------------------------------------------------------

if (verbose) message("Selecting mock samples...")
sample_metadata <- as_tibble(metadata(lupus)$sample_metadata) %>%
    mutate(across(where(is.factor), as.character))

## Keep only Healthy, European Females
samples_mock <- filter(sample_metadata,
    Sex == "Female" & Ethnicity == "European" & Status == "Healthy"
)

## Keep only the largest batches
keep_batches <- c("dmx_count_AH7TNHDMXX_YE_8-30",
                  "dmx_count_AHCM2CDMXX_YE_0831",
                  "dmx_count_BH7YT2DMXX_YE_0907")

samples_mock <- filter(samples_mock, batch_cov %in% keep_batches)

## Remove replicated samples: only 1 such case: "dmx_count_AHCM2CDMXX_YE_0831:IGTB1906_IGTB1906"
rep_sample <- samples_mock$sample[which(duplicated(samples_mock$ind_cov))]
samples_mock <- filter(samples_mock, sample != rep_sample)

## Remove Age outlier: again just 1 sample
samples_mock <- filter(samples_mock, Age < 50)

## Should have 44 patients / samples now (the two are equivalent from here on)
stopifnot(nrow(samples_mock) == 44)
stopifnot(sum(duplicated(samples_mock$ind_cov)) == 0)

lupus_mock <- lupus[, lupus$sample %in% samples_mock$sample]
colData(lupus_mock) <- droplevels(colData(lupus_mock))

## Update metadata(lupus_mock): drop patient_splits and update sample_metadata
metadata(lupus_mock)$patient_splits <- NULL
metadata(lupus_mock)$sample_metadata <- samples_mock


# Filter celltypes --------------------------------------------------------
if (verbose) message("Filtering celltypes...")

## Removing celltypes with less than 100 median cells per patient
ct_abundance <- table(lupus_mock$ind_cov, lupus_mock$ct_cov)
ct_med_abundance <- colMedians(ct_abundance)
keep_ct <- colnames(ct_abundance)[ct_med_abundance > 100]

lupus_mock <- lupus_mock[, lupus_mock$ct_cov %in% keep_ct]
colData(lupus_mock) <- droplevels(colData(lupus_mock))


# Remove non-expressed genes ----------------------------------------------
if (verbose) message("Filtering genes...")

keep_genes <- rowSums(counts(lupus_mock) > 0) >= 1
lupus_mock <- lupus_mock[keep_genes, ]

## Store mock data as new HDF5-backed SCE for faster access
if (verbose) message("Exporting data...")
lupus_mock <- saveHDF5SummarizedExperiment(
    lupus_mock, dir = h5_out_dir,
    replace = TRUE,  # to avoid errors when re-running
    verbose = TRUE
)

## counts(lupus_mock) now contains a pointer to the new (smaller) HDF5 file
stopifnot(path(counts(lupus_mock)) == file.path(h5_out_dir, "assays.h5"))

if (verbose) message("Done.")
