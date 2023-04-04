## Generates mock data replicates from the provided input data by randomly
## assigning mock labels to each patient.
## The output is a list of SCE objects for each mock replicate.

# Command line arguments --------------------------------------------------

library(argparse)
parser <- ArgumentParser()
parser$add_argument("-v", "--verbose",
    action = "store_true",
    help = "Print verbose output. Optional flag."
)
parser$add_argument("-i", "--infile",
    type = "character",
    default = "data/lupus-SCE-mock-ct-subsets/lupus-SCE-mock-T4_naive.rds",
    help = "Input data file name. Default: \'%(default)s\'."
)
parser$add_argument("-o", "--outfile",
    type = "character",
    default = "data/lupus-SCE_list-ncM-mock_replicates.rds",
    help = "Output results file name. Default: \'%(default)s\'."
)
parser$add_argument("--n_mock_replicates",
    type = "integer", default = 5L,
    help = "Integer. Number of mock data replicates. Default: %(default)s."
)

args <- parser$parse_args()
verbose <- args$verbose

if (verbose) {
    message("Using input file: ", args$infile)
    message("Using output file: ", args$outfile)
    message("Number of mock replicates: ", args$n_mock_replicates)
}

## Directory setup
here_root <- "benchmarks/lupus"
here::i_am(file.path(here_root, "scripts", "03-prepare-mock-replicates.R"))

in_file <- here::here(here_root, args$infile)
out_file <- here::here(here_root, args$outfile)

if (!file.exists(in_file)) stop("Could not find input file ", in_file)
stopifnot(dir.exists(dirname(out_file)))


## Libraries
suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(SCandwichCompanion)
})

seed <- 20220505


# Prepare data ------------------------------------------------------------

## Load data
if (verbose) message("Loading data...")
sce <- readRDS(in_file)

## Assign mock group for each patient
# Assign 50% of patients within each batch to either mock A or mock B
if (verbose) message("Creating mock groups...")
n_mock_replicates <- args$n_mock_replicates

## Repeat mock group assignment n_mock_replicates times
# TODO: these replicates only differ w.r.t. to the mock label assignment
# Hence, it is quite inefficient to replicate the entire object so many times
# In stead, we could just store a list of cell-level mock assignments.
# It does not matter too much, it just would save ~500Mb storage
set.seed(seed)
sce_list <- replicate(n_mock_replicates,
    SCandwichCompanion::createMockGroups(sce),
    simplify = FALSE
)
names(sce_list) <- paste0("replicate_", seq_len(n_mock_replicates))

## Filter genes
# uses edgeR::filterByExpr, with large.n = ncol(sce)/50 (2%) + default settings
if (verbose) message("Filtering genes...")
sce_list <- lapply(sce_list, function(x) {
    filterGenes(x, group = factor(x[["mock_group"]]))
})


# Export output -----------------------------------------------------------

if (verbose) message("Exporting output data...")
saveRDS(sce_list, file = out_file)

if (verbose) message("Done.\n")
