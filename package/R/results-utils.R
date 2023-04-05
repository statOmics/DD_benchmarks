## Get result files

#' @export
get_mock_res_files <- function(dataset,
                               methods,
                               datatype,
                               celltype = NULL,
                               n_patients = NULL) {
    out <- .get_res_files(
        dataset = dataset,
        methods = methods,
        datatype = datatype,
        which = "mock_results",
        celltype = celltype,
        n_patients = n_patients
    )
    setNames(out, methods)
}

#' @export
get_sim_res_files <- function(dataset,
                              methods,
                              datatype,
                              prop_DE,
                              celltype = NULL,
                              n_patients = NULL) {
    stopifnot(length(prop_DE) == 1)
    prop_DE_ <- sub("\\.", "_", prop_DE)
    method_prop_DE <- paste0(methods, "-prop_DE-", prop_DE_)
    out <- .get_res_files(
        dataset = dataset,
        methods = method_prop_DE,
        datatype = datatype,
        which = "sim_results",
        celltype = celltype,
        n_patients = n_patients
    )
    setNames(out, methods)
}

.get_res_files <- function(dataset,
                           methods,
                           datatype,
                           which = c("mock_results", "sim_results"),
                           celltype = NULL,
                           n_patients = NULL) {
    which <- match.arg(which)
    res_dir <- .subdir(dataset, subdir = "results")

    if (!is.null(n_patients)) dataset <- paste(dataset, n_patients, sep = "-")
    if (!is.null(celltype)) dataset <- paste(dataset, celltype, sep = "-")

    fnames <- paste(dataset, which, datatype, methods, sep = "-")
    out <- file.path(res_dir, paste0(fnames, ".rds"))
    .check_existing_files(out)

    out
}


#' @export
get_SCE_files <- function(dataset,
                          which = c("mock_replicates", "sim_replicates"),
                          celltype = NULL,
                          n_patients = NULL,
                          prop_DE = NULL) {
    which <- match.arg(which)
    stopifnot(is.null(prop_DE) || length(prop_DE) == 1)
    data_dir <- .subdir(dataset, subdir = "data")

    if (!is.null(n_patients)) dataset <- paste(dataset, n_patients, sep = "-")
    if (!is.null(celltype)) which <- paste(celltype, which, sep = "-")

    if (!is.null(prop_DE)) {
        prop_DE_ <- sub("\\.", "_", prop_DE)
        which <- paste0(which, "-prop_DE-", prop_DE_)
    }
    fnames <- paste0(dataset, "-SCE_list-", which)
    out <- file.path(data_dir, paste0(fnames, ".rds"))
    .check_existing_files(out)

    out
}


## Helper to get the requested subdirectory for a given dataset
.subdir <- function(dataset, subdir) {
    here_root <- paste0("benchmarks/", dataset)
    here::here(here_root, subdir)
}


.check_existing_files <- function(files) {
    exist <- file.exists(files)
    if (!all(exist)) {
        not_existing <- files[!exist]
        stop("The following files were not found:",
            paste0("\n  * ", not_existing),
            call. = FALSE
        )
    }
    invisible(TRUE)
}


# Extractor functions -----------------------------------------------------

## Run times for DE methods within one group (e.g. prop_DE)
#' @importFrom dplyr bind_rows
#' @export
get_runtimes <- function(x, depth) {
    out <- .extract(x, "run_time", depth = depth)
    out <- bind_rows(out, .id = "replicate")
    out <- tidyr::pivot_longer(out, everything(), names_to = "replicate")
    dplyr::rename(out, time = value)
}

## General extractor for a (nested) list of results
#' @importFrom purrr map_depth
.extract <- function(x, which, depth) {
    map_depth(x, depth, which)
}

## Get DE tables
#' @importFrom purrr map
#' @export
get_tables <- function(x, depth) {
    out <- .extract(x, "results", depth = depth)
    map(out, function(.x) {
        if (!"gene" %in% colnames(.x)) {
            .x <- tibble::rownames_to_column(.x, var = "gene")
        }
        .x
    })
}

#' @export
get_aggregated_rep_tables <- function(x, depth) {
    tables <- get_tables(x, depth = depth)
    .aggregate_rep_tables(tables)
}

## Aggregate replicate tables
.aggregate_rep_tables <- function(tables) {
    bind_rows(tables, .id = "replicate")
}

## Combine tables, using common columns
## Expects list of data.frames as input
#' @export
combine_tables <- function(tables, .id = "method", only_common_cols = TRUE) {
    all_names <- map(tables, colnames)
    if (only_common_cols) {
        common <- purrr::reduce(all_names, intersect)
        tables <- map(tables, dplyr::select, tidyselect::all_of(common))
    }
    bind_rows(tables, .id = .id)
}


# iCOBRA utils ------------------------------------------------------------

#' @export
#' @importFrom iCOBRA COBRAData calculate_adjp
prepare_COBRAData <- function(res, sce, replace_missing = TRUE) {
    ## Get p-values for iCOBRA
    pvals <- res[, c("gene", "method", "PValue")]
    pvals <- pivot_wider(pvals, names_from = method, values_from = PValue)
    pvals <- data.frame(pvals[, -1], row.names = pvals$gene)

    if (replace_missing) {
       tmp <- lapply(pvals, function(x) {
           x[is.na(x)] <- 1
           x
       })
       pvals <- data.frame(tmp, row.names = rownames(pvals))
    }
    ## Create COBRAData object with ground truth
    out <- COBRAData(
        pval = pvals,
        truth = data.frame(
            status = as.integer(rowData(sce)[["is_DE"]]),
            row.names = rownames(sce)
        )
    )
    suppressMessages(calculate_adjp(out))
}

#' @export
#' @importFrom iCOBRA fdrtpr
combine_fdrtpr_tables <- function(cobra_list) {
    out <- purrr::map_dfr(cobra_list, fdrtpr, .id = "replicate")
    dplyr::mutate(out,
           thr = as.numeric(sub("thr", "", thr)),
           method = as.factor(method)
    )
}
