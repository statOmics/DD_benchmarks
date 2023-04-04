#' Apply regular pairwise t-tests
#'
#' Perform DE analyses using pairwise t-tests as implemented in
#' [`scran::pairwiseTTests()`], blocking on batch covariates.
#'
#' @param sce A SingleCellExperiment object.
#' @param ... Ignored.
#' @param group_id Name of the `colData` column containg group identifiers.
#' @param batch_id Name of the `colData` column containg batch identifiers.
#'
#' @return Results data.frame.
#'
#' @importFrom scuttle logNormCounts
#' @importFrom S4Vectors rename
#' @export
apply_ttest <- function(sce, ...,
                        group_id = "group_id",
                        batch_id = "batch_cov") {
    sce <- logNormCounts(sce)
    grps <- colData(sce)[[group_id]]
    blocks <- factor(colData(sce)[[batch_id]])
    res <- scran::pairwiseTTests(logcounts(sce), groups = grps, block = blocks)
    out <- res$statistics[[1]]
    out$comparison <- group_id
    as.data.frame(rename(out, p.value = "PValue"))
}
