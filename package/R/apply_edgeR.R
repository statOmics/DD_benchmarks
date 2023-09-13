#' Run a GLM ...
#'
#' @param sce A SingleCellExperiment object.
#' @param coef_test Coefficient to test.
#' @param formula Design formula.
#' @param BPPARAM Parallel execution (ignored for edgeR functions)
#' @param ... Ignored.
#'
#' @return Results data.frame.
#'
#' @import edgeR
#' @export

apply_edgeR_QP <- function(sce,
                           coef_test = "group_idB",
                           formula = ~ batch_cov + group_id,
                           BPPARAM,
                           ...) {

    bin_counts <- as.matrix(assay(sce))
    cd <- colData(sce)

    y <- DGEList(counts=bin_counts, samples=cd)
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

    y <- DGEList(counts=bin_counts, samples=cd)
    design <- model.matrix(formula, y$samples)

    if(all(cd$data_type == "sce")){
        y$offset <- log(colMeans(bin_counts))
        y <- estimateDisp(y, design)
        fit <- glmQLFit(y, design, robust=T)
        res <- glmQLFTest(fit, coef = coef_test)
    } else if(all(cd$data_type == "pb")){
        y$offset <- log(cd$ncells)
        y <- estimateDisp(y, design)
        fit <- glmQLFit(y, design, robust=T)
        res <- glmQLFTest(fit, coef = coef_test)
    }
    res$table$FDR <- p.adjust(res$table$PValue, method = "BH")
    return(as.data.frame(res$table))
}

apply_edgeR_QP_optim <- function(sce,
                           coef_test = "group_idB",
                           formula = ~ batch_cov + group_id,
                           BPPARAM,
                           ...) {

    cd <- colData(sce)
    if(all(cd$data_type == "sce")){
        warning("Running apply_edgeR_QP_optim only makes sense on pseudobulk
                data: Results identical to running apply_edgeR_QP will be
                returned. ")
        res <- apply_edgeR_QP(sce,
                       coef_test = "group_idB",
                       formula = ~ batch_cov + group_id,
                       BPPARAM)
        return(res)
    }

    # Remove features detected in (nearly) all cells
    bin_counts <- as.matrix(assay(sce))
    med_detection <- rowMedians(sweep(bin_counts, 2, cd$ncells, "/"))

    y <- DGEList(counts=bin_counts[med_detection < 0.9,], samples=cd)
    design <- model.matrix(formula, y$samples)

    # Include CDR offset
    of <- colMeans(sweep(bin_counts[med_detection < 0.9,], 2, cd$ncells, "/"))
    y$offset <- log(cd$ncells*of)

    fit <- glmQLFit(y, design, robust=T, dispersion=0)

    # allow for underdispersion
    res <- glmQLFTest(fit, coef = coef_test, poisson.bound = FALSE)

    res$table$FDR <- p.adjust(res$table$PValue, method = "BH")
    return(as.data.frame(res$table))
}

apply_edgeR_NB_optim <- function(sce,
                                 coef_test = "group_idB",
                                 formula = ~ batch_cov + group_id,
                                 BPPARAM,
                                 ...) {

    cd <- colData(sce)
    if(all(cd$data_type == "sce")){
        warning("Running apply_edgeR_NB_optim only makes sense on pseudobulk
                data: Results identical to running apply_edgeR_NB will be
                returned. ")
        res <- apply_edgeR_NB(sce,
                       coef_test = "group_idB",
                       formula = ~ batch_cov + group_id,
                       BPPARAM)
        return(res)
    }

    # Remove features detected in (nearly) all cells
    bin_counts <- as.matrix(assay(sce))
    med_detection <- rowMedians(sweep(bin_counts, 2, cd$ncells, "/"))

    y <- DGEList(counts=bin_counts[med_detection < 0.9,], samples=cd)
    design <- model.matrix(formula, y$samples)

    # Include CDR offset
    of <- colMeans(sweep(bin_counts[med_detection < 0.9,], 2, cd$ncells, "/"))
    y$offset <- log(cd$ncells*of)

    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=T)

    # allow for underdispersion
    res <- glmQLFTest(fit, coef = coef_test, poisson.bound = FALSE)

    res$table$FDR <- p.adjust(res$table$PValue, method = "BH")
    return(as.data.frame(res$table))
}
