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
    default = "analysis/lupus-mock-results.Rmd",
    help = "RMarkdown file name. Default: \'%(default)s\'."
)
parser$add_argument("--out_file",
    type = "character",
    default = "analysis/lupus-mock-results.html",
    help = "Output file name. Default: \'%(default)s\'."
)
parser$add_argument("--use_celltypes",
    type = "character", nargs = "+",
    default = c("T4_naive", "B_mem", "ncM"),
    help = "Celltypes to use. Default: \'%(default)s\'"
)
parser$add_argument("--use_methods",
    type = "character", nargs = "+",
    help = "Methods to use."
)
parser$add_argument("--datatype",
    type = "character", nargs = "+",
    help = "Data type to use."
)
## Only relevant for sim-results
parser$add_argument("--prop_DE",
    type = "double", nargs = "+",
    default = NULL,
    help = "Proportion of DE to use. Default: \'%(default)s\'"
)

args <- parser$parse_args()
verbose <- args$verbose

celltypes <- args$use_celltypes
methods <- args$use_methods
datatype <- args$datatype
prop_DE <- args$prop_DE

if (verbose) {
    message("Using Rmd file: ", args$rmd_file)
    message("Using output file: ", args$out_file)
    message("Using celltypes: ", paste(celltypes, collapse = ", "))
    message("Using methods: ", paste(methods, collapse = ", "))
    message("Using datatype: ", paste(datatype, collapse = ", "))
    if (!is.null(prop_DE)) {
        message("Using prop_DE: ", paste(prop_DE, collapse = ", "))
    }
}


## Directory setup
here_root <- "benchmarks/lupus"
here::i_am(file.path(here_root, "scripts", "05-render-lupus-results.R"))

rmd_file <- here::here(here_root, args$rmd_file)
out_file <- here::here(here_root, args$out_file)

if (verbose) {
    message("Rendering output to: ", out_file)
}

stopifnot(file.exists(rmd_file))


# Render Rmd --------------------------------------------------------------

rmd_params <- list(methods = methods,
                   celltypes = celltypes,
                   datatype = datatype)
if (!is.null(prop_DE)) {
    rmd_params[["prop_DE"]] <- prop_DE
}

rmarkdown::render(rmd_file,
    output_file = out_file,
    params = rmd_params
)
