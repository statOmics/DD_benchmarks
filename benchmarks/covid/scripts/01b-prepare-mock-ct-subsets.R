## Create single-celltype subsets from the mock data
## =================================================

# Command line arguments --------------------------------------------------

library(argparse)
parser <- ArgumentParser()
parser$add_argument("-v", "--verbose",
    action = "store_true",
    help = "Print verbose output. Optional flag."
)
parser$add_argument("-i", "--infile",
    type = "character",
    default = "data/covid-SCE-prepared.rds",
    help = "Input file name for the SummarizedExperiment. Default: \'%(default)s\'."
)
parser$add_argument("--out_dir",
    type = "character",
    default = "data/covid-SCE-mock-ct-subsets",
    help = "Output directory name for the generated data. Default: \'%(default)s\'."
)
parser$add_argument("--use_celltypes",
    type = "character", nargs = "+",
    default = c("class_switched_memory_B_cell_Critical",
                "immature_B_cell_Healthy",
                "naive_B_cell_Healthy"),
    help = "Celltypes to use. Default: \'%(default)s\'"
)

args <- parser$parse_args()
verbose <- args$verbose

use_celltypes <- args$use_celltypes

if (verbose) {
    message("Using input data: ", args$infile)
    message("Using output directory: ", args$out_dir)
    message("Using celltypes: ", paste(use_celltypes, collapse = ", "))
}

## Directory setup
here_root <- "benchmarks/covid"
here::i_am(file.path(here_root, "scripts", "01b-prepare-mock-ct-subsets.R"))

in_file <- here::here(here_root, args$infile)
out_dir <- here::here(here_root, args$out_dir)

if (!file.exists(in_file)) stop("Could not find file ", in_file)
stopifnot(dir.exists(dirname(out_dir)))

## Libraries
suppressPackageStartupMessages({
    library(SingleCellExperiment)
})

if (verbose) message("Loading data...")
covid_mock <- readRDS(file = in_file)

# Make individual celltype subsets ----------------------------------------
if (verbose) message("Making subsets...")

missing <- !(use_celltypes %in% covid_mock$ct_cov)
if (any(missing)) {
    stop("Celltypes ", paste(use_celltypes[missing], collapse = ", "), " not found.")
}

## Creating separate SCE subsets with in-memory counts for easier later access
ct_subs <- lapply(use_celltypes, function(x) {
    sce <- covid_mock[, covid_mock$ct_cov == x]
    colData(sce) <- droplevels(colData(sce))
    counts(sce) <- as(counts(sce), "sparseMatrix")
    keep_genes <- rowSums(counts(sce) > 0) >= 1
    sce <- sce[keep_genes, ]
    sce
})
names(ct_subs) <- use_celltypes

## Export subsets
if (verbose) message("Exporting data...")

if (dir.exists(out_dir)) fs::dir_delete(out_dir) # clean up during re-runs
fs::dir_create(out_dir) # ensures timestamp is updated

lapply(names(ct_subs), function(ct_name) {
    x <- ct_subs[[ct_name]]
    out_file <- file.path(out_dir, paste0("covid-SCE-mock-", ct_name, ".rds"))
    saveRDS(x, out_file)
    out_file
})

if (verbose) message("Done.")
