# Command line arguments --------------------------------------------------

library(argparse)
parser <- ArgumentParser()
parser$add_argument("-v", "--verbose",
    action = "store_true",
    help = "Print verbose output. Optional flag."
)
parser$add_argument("-i", "--infile",
    type = "character",
    default = "../lupus/data/lupus-SCE_list-ncM-mock_replicates.rds",
    help = "Input data file name. Default: \'%(default)s\'."
)
parser$add_argument("-o", "--outfile",
    type = "character",
    default = "data/lupus-n_patients-10-SCE_list-ncM-mock_replicates.rds",
    help = "Output results file name. Default: \'%(default)s\'."
)
parser$add_argument("--n_patients",
    type = "integer", default = 10,
    help = "Integer. Number of patients to use. If missing, use all patients."
)

args <- parser$parse_args()
verbose <- args$verbose

if (verbose) {
    message("Using input file: ", args$infile)
    message("Using output file: ", args$outfile)
    message("Number of patients: ", args$n_patients)
}

## Directory setup
here::i_am(file.path("benchmarks", "lupus-n_patients/scripts", "01-sample-patients.R"))

in_file <- sub("^\\.\\./", "", args$infile) # Resolve relative paths
in_file <- here::here("benchmarks", in_file)
out_file <- here::here("benchmarks/lupus-n_patients", args$outfile)

if (!file.exists(in_file)) stop("Could not find input file ", in_file)
stopifnot(dir.exists(dirname(out_file)))

n_patients <- args$n_patients

## Libraries
suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(DDCompanion)
    library(scuttle)
})

## Load data
if (verbose) message("Loading data...")
sce_list <- readRDS(in_file)

## Check for expected data structure
stopifnot(all(unlist(lapply(sce_list, is, "SingleCellExperiment"))))

seed <- 20220426


# Helpers -----------------------------------------------------------------

samplePatients <- function(x, n_patients) {
    patients <- as.character(colData(x)[["ind_cov"]])
    batches <- as.character(colData(x)[["batch_cov"]])

    patients_by_batch <- lapply(split(patients, batches), unique)

    ## Only use as many batches as needed to reach the desired n_patients
    tmp <- sort(lengths(patients_by_batch), decreasing = TRUE)
    use_batches <- names(tmp)[1:which(n_patients <= cumsum(tmp))[1]]
    x <- x[, x$batch_cov %in% use_batches]
    colData(x) <- droplevels(colData(x))

    ## No need to sample if we already have the desired number of patients
    if (sum(lengths(patients_by_batch[use_batches])) > n_patients) {
        ## Stratified sampling of subjects
        patients <- as.character(colData(x)[["ind_cov"]])
        mock_groups <- as.character(colData(x)[["mock_group"]])

        patients_by_group <- lapply(split(patients, mock_groups), unique)

        n_groups <- length(patients_by_group)
        keep_patients <- lapply(patients_by_group, sample, size = n_patients / n_groups)
        x <- x[, x$ind_cov %in% unlist(keep_patients)]
        colData(x) <- droplevels(colData(x))
    }

    ## Remove unexpressed genes and return
    x[rowSums(counts(x) > 0) >= 1, ]
}


# Sub-sample patients -----------------------------------------------------

if (verbose) message("Subsampling patients...")
set.seed(seed)
sce_list <- lapply(sce_list, samplePatients, n_patients = n_patients)

## Filter genes
if (verbose) message("Filtering genes...")
sce_list <- lapply(sce_list, function(x) {
    filterGenes(x, group = factor(x[["mock_group"]]))
})

## Binarize the data
sce_list_bin <- lapply(sce_list, function(element){
    assay(element)[assay(element)>=1] <- 1
    element$data_type <- "sce"
    return(element)
})

## Pseudobulk the data
pb_list_bin <- lapply(sce_list_bin, function(element){
    pb_bin <- aggregateAcrossCells(element,
                                   ids = element$ind_cov)
    colnames(pb_bin) <- paste0("identifier_", 1:ncol(pb_bin))
    pb_bin$data_type <- "pb"
    return(pb_bin)
})

# Export output -----------------------------------------------------------

if (verbose) message("Exporting output data...")
saveRDS(pb_list_bin, file = out_file)

if (verbose) message("Done.\n")
