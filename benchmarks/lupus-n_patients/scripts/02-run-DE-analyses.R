# Command line arguments --------------------------------------------------

library(argparse)
parser <- ArgumentParser()
parser$add_argument("-v", "--verbose",
    action = "store_true",
    help = "Print verbose output. Optional flag."
)
parser$add_argument("-i", "--infile",
    type = "character",
    default = "data/lupus-n_patients-10-SCE_list-ncM-mock_replicates.rds",
    help = "Input data file name. Default: \'%(default)s\'."
)
parser$add_argument("-o", "--outfile",
    type = "character",
    default = "results/lupus-n_patients-10-ncM-mock_results-muscat.rds",
    help = "Output results file name. Default: \'%(default)s\'."
)
parser$add_argument("--method",
    type = "character", default = "muscat",
    help = "DE methods to be run on the input data. Default: \'%(default)s\'."
)

args <- parser$parse_args()
verbose <- args$verbose

if (verbose) {
    message("Using input file: ", args$infile)
    message("Using output file: ", args$outfile)
    message("Using DE method: ", args$method)
}

## Directory setup
here_root <- "benchmarks/lupus-n_patients"
here::i_am(file.path(here_root, "scripts", "02-run-DE-analyses.R"))

in_file <- here::here(here_root, args$infile)
out_file <- here::here(here_root, args$outfile)

if (!file.exists(in_file)) stop("Could not find input file ", in_file)
stopifnot(dir.exists(dirname(out_file)))


## Libraries
suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(purrr)
    library(DDCompanion)
})


# Prepare data ------------------------------------------------------------

## Load data
if (verbose) message("Loading data...")
pb_list_bin <- readRDS(in_file)

## Check for expected data structure
stopifnot(all(unlist(map(pb_list_bin, is, "SingleCellExperiment"))))

pb_list_bin <- map(pb_list_bin, prepare_SCE,
    group_id = "mock_group",
    subject_id = "ind_cov",
    cluster_id = "ct_cov"
)


# Run DE analyses ---------------------------------------------------------

## Arguments for DDCompanion::run_de_method()
tmp <- pb_list_bin[[1]] # same for all replicates
if (nlevels(tmp$batch_cov) == 1) {
    formula <- ~group_id
} else {
    formula <- ~ batch_cov + group_id
}

## Loop over list and run DE method on each SCE
if (verbose) message("Runnig DE analyses...")
out <- map(pb_list_bin, run_de_method,
    method = args$method,
    formula = formula,
    coef = "group_idB"
)

## Export results
if (verbose) message("Exporting results...")
saveRDS(out, file = out_file)

if (verbose) message("Done.\n")
