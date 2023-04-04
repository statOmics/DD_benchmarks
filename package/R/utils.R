#' Create subject-level mock groups
#'
#' Randomly assign subjects to mock groups for a single-cell experiment.
#'
#' @param x A \linkS4class{SummarizedExperiment} object.
#' @param labs Mock labels to use. The length determines the number of groups.
#' @param subject_id `colData` column name that contains the subject labels.
#' @param batch_id `colData` column name that contains batch labels.
#'
#' @return
#' An object of class `class(x)` with mock labels added in the `"mock_group"`
#' column of the `colData`.
#' @export
createMockGroups <- function(x, labs = c("A", "B"),
                             subject_id = "ind_cov",
                             batch_id = "batch_cov") {

    ## Avoid factor shenanigans when assigning
    subjects <- as.character(colData(x)[[subject_id]])
    batches <- as.character(colData(x)[[batch_id]])

    ## Sample mock groups stratified by batch
    batch_subjects <- split(subjects, batches)
    ## Need to remove batch names to avoid problems later with `unlist()`
    batch_subjects <- unname(lapply(batch_subjects, unique))

    mock_labs <- lapply(batch_subjects, function(subjects) {
        ## Ensures we have balanced groups
        out <- rep(labs, each = ceiling(length(subjects) / length(labs)))
        out <- sample(out, size = length(subjects)) # randomly permute
        names(out) <- subjects
        out
    })
    mock_labs <- unlist(mock_labs)

    colData(x)[["mock_group"]] <- mock_labs[subjects]
    x
}


#' Filter genes for mock analyses
#'
#' Uses `edgeR::filterByExpr()` to remove lowly expressed genes, taking the
#' design of the experiment into account.
#'
#' @param sce A \linkS4class{SingleCellExperiment} object.
#' @inheritParams edgeR::filterByExpr
#' @inheritDotParams edgeR::filterByExpr
#'
#' @return A \linkS4class{SingleCellExperiment} object with filtered genes removed.
#' @export
filterGenes <- function(sce, group,
                        min.count = 1L, min.total.count = 10L,
                        ...) {
    large_n <- round(ncol(sce) / 50)  # 2% of cells
    keep_genes <- edgeR::filterByExpr(
        y = counts(sce), group = group,
        min.count = min.count, min.total.count = min.total.count,
        large.n = large_n, min.prop = 0,
        ...
    )
    sce[keep_genes, ]
}


#' Sub-sample cells per subject
#'
#' Down-samples the number of cells for each subject to a given number. Then
#' calls \code{\link{filterGenes}} on each subset.
#'
#' @param x A \linkS4class{SummarizedExperiment} object.
#' @param n_cells Number of cells per subject to keep. Possible to supply more
#'   than one value.
#' @param subject_id `colData` column name containing the subject labels
#' @param group_id `colData` column name containing the group labels. Used for
#'   gene filtering with `filterGenes()`.
#'
#' @return A list of objects of class `class(x)` with length `length(n_cells)`.
#'   Each cobject coresponds to a down-sampled dataset
#' @export
sampleCellsPerSubject <- function(x, n_cells,
                                  subject_id = "ind_cov",
                                  group_id = "mock_group") {
    cd <- colData(x)
    split_idx <- split(seq_len(nrow(cd)), factor(cd[[subject_id]]))

    out <- vector("list", length = length(n_cells))

    for (i in seq_along(n_cells)) {
        n <- n_cells[[i]]
        keep_cells <- lapply(split_idx, function(y) {
            if (length(y) <= n) {
                y
            } else {
                sample(y, size = n, replace = FALSE)
            }
        })
        x <- x[, unlist(keep_cells)]
        ## Remove unexpressed genes
        out[[i]] <- filterGenes(x, group = factor(cd[[group_id]]))
    }
    names(out) <- paste0("cells_per_subject_", n_cells)
    out
}
