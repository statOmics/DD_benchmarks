## Render lupus-results.Rmd with input parameters
## ==============================================


# Command line arguments --------------------------------------------------

library(argparse)
parser <- ArgumentParser()
parser$add_argument("-v", "--verbose",
    action = "store_true",
    help = "Print verbose output. Optional flag."
)
parser$add_argument("--rmd-file",
    type = "character",
    default = "analysis/lupus-n_patients-mock-results.Rmd",
    help = "RMarkdown file name. Default: \'%(default)s\'."
)
parser$add_argument("--out_file",
    type = "character",
    default = "analysis/lupus-n_patients-mock-results-ncM.html",
    help = "Output file name. Default: \'%(default)s\'."
)
parser$add_argument("--n_patients",
    type = "integer", nargs = "+", # allow 0 or 1 argument
    default = 10,
    help = "Integer. Number of patients to use. If missing, use all patients."
)
parser$add_argument("--use_celltype",
    type = "character",
    default = "ncM",
    help = "Celltype to use. Default: \'%(default)s\'"
)
parser$add_argument("--use_methods",
    type = "character", nargs = "+",
    default = c("muscat"),
    help = "Methods to use."
)

args <- parser$parse_args()
verbose <- args$verbose

n_patients <- args$n_patients

celltype <- args$use_celltype
methods <- args$use_methods

if (verbose) {
    message("Using Rmd file: ", args$rmd_file)
    message("Using output file: ", args$out_file)
    message("Number of patients: ", args$n_patients)
    message("Using celltype: ", paste(celltype, collapse = ", "))
    message("Using methods: ", paste(methods, collapse = ", "))
}


## Directory setup
here_root <- "benchmarks/lupus-n_patients"
here::i_am(file.path(here_root, "scripts", "03-render-lupus-n_patients-results.R"))

rmd_file <- here::here(here_root, args$rmd_file)
out_file <- here::here(here_root, args$out_file)

if (verbose) {
    message("Rendering output to: ", out_file)
}

stopifnot(file.exists(rmd_file))


# Render Rmd --------------------------------------------------------------

rmarkdown::render(rmd_file,
    output_file = out_file,
    params = list(
        n_patients = n_patients,
        methods = methods,
        celltype = celltype
    )
)
