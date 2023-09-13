#' @export
run_de_method <- function(sce, method, formula, coef, BPPARAM, combined = FALSE) {

    de_fun <- DDCompanion:::.get_method(method)

    ## Run DE method and measure elapsed runtime
    t <- system.time({
        res <- de_fun(sce, coef = coef, formula = formula, BPPARAM = BPPARAM)
    })[[3]]
    out <- list(results = res, run_time = t)
    out
}

#' @export
prepare_SCE <- function(sce, group_id, subject_id, cluster_id) {
    ## Ensure consistent covariate names
    stopifnot(!is.null(colData(sce)[[group_id]]))
    stopifnot(!is.null(colData(sce)[[subject_id]]))
    stopifnot(!is.null(colData(sce)[[cluster_id]]))

    sce[["group_id"]] <- sce[[group_id]]
    sce[["subject_id"]] <- sce[[subject_id]]
    sce[["cluster_id"]] <- sce[[cluster_id]]

    sce
}

.get_method <- function(method) {
    utils::getFromNamespace(paste0("apply_", method), "DDCompanion")
}
