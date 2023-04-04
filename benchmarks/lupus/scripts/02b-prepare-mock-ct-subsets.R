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
    default = "data/lupus-SCE-mock-h5/se.rds",
    help = "Input file name for the SummarizedExperiment. Default: \'%(default)s\'."
)
parser$add_argument("--out_dir",
    type = "character",
    default = "data/lupus-SCE-mock-ct-subsets",
    help = "Output directory name for the H5 files. Default: \'%(default)s\'."
)
parser$add_argument("--use_celltypes",
    type = "character", nargs = "+",
    default = c("T4_naive", "B_mem", "ncM"),
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
here_root <- "benchmarks/lupus"
here::i_am(file.path(here_root, "scripts", "02b-prepare-mock-ct-subsets.R"))

h5_in_dir <- here::here(here_root, dirname(args$infile))
out_dir <- here::here(here_root, args$out_dir)

if (!dir.exists(h5_in_dir)) stop("Could not find input directory ", h5_in_dir)
stopifnot(dir.exists(dirname(out_dir)))

## Libraries
suppressPackageStartupMessages({
    library(HDF5Array)
    library(SingleCellExperiment)
})

if (verbose) message("Loading data...")
lupus_mock <- loadHDF5SummarizedExperiment(dir = h5_in_dir)


# Make individual celltype subsets ----------------------------------------
if (verbose) message("Making subsets...")

missing <- !(use_celltypes %in% lupus_mock$ct_cov)
if (any(missing)) {
    stop("Celltypes ", paste(use_celltypes[missing], collapse = ", "), " not found.")
}

## Creating separate SCE subsets with in-memory counts for easier later access
ct_subs <- lapply(use_celltypes, function(x) {
    sce <- lupus_mock[, lupus_mock$ct_cov == x]
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
    out_file <- file.path(out_dir, paste0("lupus-SCE-mock-", ct_name, ".rds"))
    saveRDS(x, out_file)
    out_file
})

if (verbose) message("Done.")
