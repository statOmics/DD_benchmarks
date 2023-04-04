#' @export
run_de_method <- function(sce, method, formula, coef, combined = FALSE) {

    de_fun <- .get_method(method)

    if (combined) {
       out <- .run_de_method_combined(
           sce = sce, de_fun = de_fun,
           formula = formula, coef = coef
       )
    } else {
        ## Run DE method and measure elapsed runtime
        t <- system.time({
            res <- de_fun(sce, coef = coef, formula = formula)
        })[[3]]
        out <- list(results = res, run_time = t)
    }

    out
}

#' @importFrom stats relevel
.run_de_method_combined <- function(sce, de_fun, formula, coef) {
    ## Test the mock group effect within each celltype
    ## We'll test the main group_idB effect and just relevel so that each
    ## celltype is taken as the reference once
    clusters <- levels(sce$cluster_id)

    ## Run DE method and measure elapsed runtime
    out <- vector("list", length = length(clusters))
    names(out) <- clusters
    for (ct in clusters) {
        ## Relevel
        sce$cluster_id <- relevel(sce$cluster_id, ref = ct)
        t <- system.time({
            res <- de_fun(sce,
                coef = coef,
                formula = formula, combined = TRUE
            )
        })[[3]]
        out[[ct]] <- list(results = res, run_time = t)
    }
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
    utils::getFromNamespace(paste0("apply_", method), "SCandwichCompanion")
}
