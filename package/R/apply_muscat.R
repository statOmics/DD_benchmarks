#' Apply muscat with default parameters
#'
#' When analysing multiple celltypes combined, uses a workaround with
#' `scran::pseudobulkDGE()` to apply the pseudobulk DE approach.
#' `edgeR::glmQLFTest()` is used as the underlying DE method.
#'
#' @param sce A SingleCellExperiment object.
#' @param coef Coefficient to test.
#' @param formula Design formula.
#' @param subject_id Name of the `colData` column containg subject identifiers.
#' @param group_id Name of the `colData` column containg group identifiers.
#' @param cluster_id Name of the `colData` column containg cluster identifiers.
#' @param ... Ignored.
#' @param combined Logical. Should the pseudobulk approach be applied on
#'     multiple celltypes together?. Default: `FALSE`.
#'
#' @return Results data.frame.
#'
#' @importFrom muscat prepSCE
#' @importFrom stats model.matrix
#' @importFrom S4Vectors rename
#' @export
apply_muscat <- function(sce, coef = "group_idB",
                         formula = ~ batch_cov + group_id,
                         subject_id = "subject_id",
                         group_id = "group_id",
                         cluster_id = "cluster_id", ...,
                         combined = FALSE) {

    ## muscat prep
    sce <- prepSCE(sce, gid = group_id, sid = subject_id, kid = cluster_id)

    if (combined) {
        .apply_pseudoBulkDGE(sce = sce, coef, formula = formula)
    } else {
        .apply_muscat(sce = sce, coef, formula = formula, filter = "none")
    }
}


## muscat helper
#' @importFrom muscat aggregateData pbDS
.apply_muscat <- function(sce, coef, formula, filter = "both") {
    pb <- aggregateData(sce,
        assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"),
        scale = FALSE, verbose = FALSE
    )

    ## Make design and make coef numeric
    design <- model.matrix(formula, data = colData(pb))
    coef <- match(coef, colnames(design))
    if (is.na(coef)) stop("`coef = `", coef, " not found in design.")
    ## Run pbDS with default parameters
    res <- pbDS(pb,
        method = "edgeR", design = design, coef = coef,
        min_cells = 10, filter = filter, verbose = FALSE
    )
    ## Expected to be run on single celltype and single coef
    out <- res$table[[1]][[1]]
    out <- rename(out,
        p_val = "PValue",
        p_adj.loc = "FDR",
        coef = "comparison"
    )
    out$p_adj.glb <- NULL # to avoid confusion

    out
}

## Use scran::pseudoBulkDGE() to apply pseudobulk (muscat) approach to combined
## celltypes
#' @importFrom scran pseudoBulkDGE
.apply_pseudoBulkDGE <- function(sce, coef, formula) {
    ## Aggregate data per sample-cluster combination
    ids <- colData(sce)[c("sample_id", "cluster_id")]
    pb <- scuttle::aggregateAcrossCells(sce, ids = ids)
    out <- pseudoBulkDGE(pb,
        label = rep("all", ncol(pb)),
        design = formula, coef = coef,
        method = "edgeR", robust = FALSE, # for consistency with muscat
        condition = pb$group_id
    )
    out <- out$all
    out$gene <- rownames(out)
    out$comparison <- coef
    as.data.frame(out)
}
