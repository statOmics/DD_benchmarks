#' Run a GLM ...
#'
#' @param sce A SingleCellExperiment object.
#' @param coef_test Coefficient to test.
#' @param formula Design formula.
#' @param ... Ignored.
#'
#' @return Results data.frame.
#'
#' @import edgeR
#' @export

apply_edgeR_QP <- function(sce,
                           coef_test = "group_idB",
                           formula = ~ batch_cov + group_id,
                           ...) {

    bin_counts <- as.matrix(assay(sce))
    cd <- colData(sce)

    y <- DGEList(counts=bin_counts[1:100,], samples=cd)
    design <- model.matrix(formula, y$samples)

    if(all(cd$data_type == "sce")){
        y$offset <- log(colMeans(bin_counts))
        fit <- glmQLFit(y, design, robust=T, dispersion=0)
        res <- glmQLFTest(fit, coef = coef_test)
    } else if(all(cd$data_type == "pb")){
        y$offset <- log(cd$ncells)
        fit <- glmQLFit(y, design, robust=T, dispersion=0)
        res <- glmQLFTest(fit, coef = coef_test)
    }
    res$table$FDR <- p.adjust(res$table$PValue, method = "BH")
    return(as.data.frame(res$table))
}

apply_edgeR_NB <- function(sce,
                           coef_test = "group_idB",
                           formula = ~ batch_cov + group_id,
                           ...) {

    bin_counts <- as.matrix(assay(sce))
    cd <- colData(sce)

    y <- DGEList(counts=bin_counts[1:100,], samples=cd)
    design <- model.matrix(formula, y$samples)

    if(all(cd$data_type == "sce")){
        y$offset <- log(colMeans(bin_counts))
        y <- estimateDisp(y, design)
        fit <- glmQLFit(y, design, robust=T, dispersion=0)
        res <- glmQLFTest(fit, coef = coef_test)
    } else if(all(cd$data_type == "pb")){
        y$offset <- log(cd$ncells)
        y <- estimateDisp(y, design)
        fit <- glmQLFit(y, design, robust=T, dispersion=0)
        res <- glmQLFTest(fit, coef = coef_test)
    }
    res$table$FDR <- p.adjust(res$table$PValue, method = "BH")
    return(as.data.frame(res$table))
}
