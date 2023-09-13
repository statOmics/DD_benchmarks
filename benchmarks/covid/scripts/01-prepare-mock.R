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
    default = "data/covid-SCE-cleaned.rds",
    help = "Input data file name. Default: \'%(default)s\'."
)
parser$add_argument("-o", "--outfile",
    type = "character",
    default = "data/covid-SCE-prepared.rds",
    help = "Output file name for the SummarizedExperiment. Default: \'%(default)s\'."
)

args <- parser$parse_args()
verbose <- args$verbose

if (verbose) {
    message("Using input file: ", args$infile)
    message("Using output file: ", args$outfile)
}

## Directory setup
here_root <- "benchmarks/covid"
here::i_am(file.path(here_root, "scripts", "01-prepare-mock.R"))

in_file <- here::here(here_root, args$infile)
out_file <- here::here(here_root, args$outfile)

if (!file.exists(in_file)) stop("Could not find input file ", in_file)
stopifnot(dir.exists(dirname(out_file)))

## Libraries
suppressPackageStartupMessages({
    library(scuttle)
    library(dplyr)
})

if (verbose) message("Loading data...")
covid <- readRDS(in_file)

seed <- 20230526

# Select samples ----------------------------------------------------------

if (verbose) message("Selecting samples for mock analysis...")
cell_metadata <- colData(covid)

# Wrangling of cell-level data
covid$sex <- as.factor(covid$sex)
covid$cell_type_curated <- gsub(" ", "_", covid$cell_type_curated)
covid$cell_type_curated <- as.factor(covid$cell_type_curated)
covid$site_sex <- as.factor(paste0(covid$Site, "_", covid$sex)) # batch effect

# Remove patients from the Sanger sequencing site
covid <- covid[,covid$Site != "Sanger"]

# Only retain class switched memory B cells, immature B cells and naive B cells
covid <- covid[,covid$cell_type_curated %in% c("class_switched_memory_B_cell",
                                               "immature_B_cell",
                                               "naive_B_cell")]

# Only retain class healthy, mildly diseased, moderately diseased and critically
# diseased patients. For these benchmarks, there are too few patients in the
# asymptomatic and severely diseased patients groups to properly test what we
# want to test.
covid <- covid[,covid$Status_on_day_collection_summary %in% c("Healthy",
                                                              "Mild",
                                                              "Moderate",
                                                              "Critical")]

# Remove patient-celltype combinations with less than 20 cells (as per the
# original publication for this dataset DOI: 10.1126/science.abf197)
covid$filtersample <- paste0(covid$donor_id, "_", covid$cell_type_curated)
covid <- covid[, covid$filtersample %in% names(which(table(covid$filtersample) >= 20))]

# I will treat cell-type/status combinations as cell types. I.e., this will
# allow me to easily perform mock analyses and simulations for each cell type,
# but separately for each status.
covid$ct_cov <- paste0(covid$cell_type_curated, "_", covid$Status_on_day_collection_summary)

# Remove non-expressed genes ----------------------------------------------
if (verbose) message("Filtering genes...")

keep_genes <- rowSums(counts(covid) > 0) >= 1
covid <- covid[keep_genes, ]

## Store mock data
if (verbose) message("Exporting data...")
saveRDS(covid, out_file)

if (verbose) message("Done.")
