#' Simulate DE by swapping genes
#'
#' Simple DE simulator that takes a dataset with a grouping factor and induces
#' DE by swapping some genes in one of the groups but not the other. This should
#' ensure keeping the same data structure as the original data without having to
#' estimate distribution parameters or making assumptions.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param groups A vector of length equal to `ncol(x)` specifying the group to
#'   which each cell is assigned. DE will be induced between these groups.
#' @param prop_DE Numeric scalar specifying the proportion of genes that will be
#'   simulated as DE. Default: 1%.
#'
#' @return
#' A \linkS4class{SingleCellExperiment} object with DE induced between the
#' specified `groups`. The simulated counts are contained within the `"counts"`
#' assay of the returned object. The original `colData(x)` and `rowData(x)` are
#' retained with additional columns in the latter:
#'
#' * `"is_DE"`: logical vector indicating the ground truth status of each gene
#' * `"swapped_gene"`: character vector indicating the original gene with which
#' the gene was swapped
#'
#' Note that normalized values in `logcounts(x)` are not retained. For simple
#' library size-factor normalization, size factors for each cell will be
#' identical before and after simulation. However, for more sophisticated
#' methods, such as normalization by deconvolution, the size factors might
#' differ due to the gene swapping.
#'
#' @examples
#' set.seed(123)
#' sce <- mockMultiSubjectSCE()
#'
#' sim_sce <- simulate_DE(sce, groups = sce$mock_group, prop_DE = 0.1)
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom methods is
#' @export
simulate_DE <- function(x, groups, prop_DE = 0.01) {
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(length(groups) == ncol(x))

    sel_group <- sample(unique(groups), 1)
    sel_group_idx <- which(groups == sel_group)

    n_rows <- nrow(x)
    n_DE <- round(prop_DE * n_rows)
    stopifnot(n_DE > 1L)
    de_genes <- sample(n_rows, size = n_DE)

    ## Scramble genes only within the selected group
    scrambling <- .scramble_rows(de_genes)
    cnts <- counts(x)
    sim_cnts <- cnts
    sim_cnts[de_genes, sel_group_idx] <- cnts[scrambling$x, sel_group_idx]

    ## It's possible that a single gene remains unswapped
    row_data <- cbind(rowData(x), is_DE = logical(n_rows))
    row_data[de_genes, "is_DE"] <- scrambling$swapped
    row_data[["swapped_gene"]] <- rep(NA_character_, n_rows)
    row_data[de_genes, "swapped_gene"] <- rownames(x)[scrambling$x]

    SingleCellExperiment(
        assays = list(counts = sim_cnts),
        colData = colData(x),
        rowData = row_data
    )
}

## Helper to scramble rows of a matrix
.scramble_rows <- function(rows) {
    ## Randomly permute rows
    scrambled <- sample(rows)
    non_swapped <- rows == scrambled

    ## Keep scrambling until all but 1 are
    while (sum(non_swapped) > 1L) {
        if (sum(non_swapped) == 2L) {
            ## Just swap if there are 2 left
            scrambled[non_swapped] <- rev(scrambled[non_swapped])
        } else {
            scrambled[non_swapped] <- sample(scrambled[non_swapped])
        }
        non_swapped <- rows == scrambled
    }

    list(x = scrambled, swapped = !non_swapped)
}
