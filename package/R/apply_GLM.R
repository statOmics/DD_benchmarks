#' Run a GLM ...
#'
#' @param sce A SingleCellExperiment object.
#' @param coef_test Coefficient to test.
#' @param formula Design formula.
#' @param ... Ignored.
#'
#' @return Results data.frame.
#'
#' @importFrom limma squeezeVar
#' @export

apply_bGLM <- function(sce,
                       coef_test = "group_idB",
                       formula = ~ batch_cov + group_id,
                       ...) {

    # Prepare for fit
    bin_counts <- as.matrix(assay(sce))
    cd <- colData(sce)

    if(all(cd$data_type == "sce")){

        formula <- as.formula(paste("y", paste(formula, collapse = " ")))

        # fit models
        mod <- lapply(seq_len(nrow(bin_counts))[1:100], function(i){
            y <- bin_counts[i,]
            cd$y <- y
            mod <- glm(formula = formula, family = "binomial", data = cd)
            return(mod)
        })

    } else if(all(cd$data_type == "pb")){

        full_formula <- as.formula(paste("cbind(succes,failure)",
                                         paste(formula, collapse = " ")))

        mod <- lapply(seq_len(nrow(bin_counts))[1:100], function(i){
            cd$succes <- bin_counts[i,]
            cd$failure <- cd$ncells - cd$succes
            mod <- glm(formula = full_formula,
                       family = "binomial",
                       data = cd)
            return(mod)
        })
    }

    # test coefficient of interest
    sum_mod <- lapply(mod, summary)
    beta <- sapply(sum_mod, function(summod){coef(summod)[coef_test,1]})
    var_unscaled <- sapply(sum_mod, function(summod){summod$cov.unscaled[coef_test,2]})
    SE <- sapply(sum_mod, function(summod){coef(summod)[coef_test,2]})
    disp <- sapply(sum_mod, function(summod){summod$dispersion})
    pval <- sapply(sum_mod, function(summod){coef(summod)[coef_test,4]})
    res <- data.frame(gene = rownames(bin_counts)[1:100],
                      cluster_id = levels(cd$cluster_id),
                      logOdds = beta,
                      var_unscaled = var_unscaled,
                      dispersion = disp,
                      SE = SE,
                      Wald = (beta^2)/SE,
                      PValue = pval,
                      FDR = p.adjust(pval, method = "BH"),
                      comparsion = coef_test)
    return(res)
}


apply_qbGLM <- function(sce,
                        coef_test = "group_idB",
                        formula = ~ batch_cov + group_id,
                        ...) {

    # Prepare for fit
    bin_counts <- as.matrix(assay(sce))
    cd <- colData(sce)

    if(all(cd$data_type == "sce")){

        formula <- as.formula(paste("y", paste(formula, collapse = " ")))

        # fit models
        mod <- lapply(seq_len(nrow(bin_counts))[1:100], function(i){
            y <- bin_counts[i,]
            cd$y <- y
            mod <- glm(formula = formula, family = "quasibinomial", data = cd)
            return(mod)
        })

    } else if(all(cd$data_type == "pb")){

        full_formula <- as.formula(paste("cbind(succes,failure)",
                                         paste(formula, collapse = " ")))

        mod <- lapply(seq_len(nrow(bin_counts))[1:100], function(i){
            cd$succes <- bin_counts[i,]
            cd$failure <- cd$ncells - cd$succes
            mod <- glm(formula = full_formula,
                       family = "quasibinomial",
                       data = cd)
            return(mod)
        })
    }

    # test coefficient of interest
    sum_mod <- lapply(mod, summary)
    beta <- sapply(sum_mod, function(summod){coef(summod)[coef_test,1]})
    var_unscaled <- sapply(sum_mod, function(summod){summod$cov.unscaled[coef_test,2]})
    SE <- sapply(sum_mod, function(summod){coef(summod)[coef_test,2]})
    disp <- sapply(sum_mod, function(summod){summod$dispersion})
    pval <- sapply(sum_mod, function(summod){coef(summod)[coef_test,4]})
    res <- data.frame(gene = rownames(bin_counts)[1:100],
                      cluster_id = levels(cd$cluster_id),
                      logOdds = beta,
                      var_unscaled = var_unscaled,
                      dispersion = disp,
                      SE = SE,
                      Wald = (beta^2)/SE,
                      PValue = pval,
                      FDR = p.adjust(pval, method = "BH"),
                      comparsion = coef_test)
    return(res)
}

apply_qbGLM_offset <- function(sce,
                               coef_test = "group_idB",
                               formula = ~ batch_cov + group_id,
                               ...) {

    # Prepare for fit
    bin_counts <- as.matrix(assay(sce))
    cd <- colData(sce)

    if(all(cd$data_type == "sce")){
        offset <- colMeans(bin_counts)
        cd$logitOffset <- log(offset/(1-offset))

        formula <- as.formula(paste("y", paste(formula, collapse = " "),
                                    "+ offset(logitOffset)"))

        # fit models
        mod <- lapply(seq_len(nrow(bin_counts))[1:100], function(i){
            y <- bin_counts[i,]
            cd$y <- y
            mod <- glm(formula = formula, family = "quasibinomial", data = cd)
            return(mod)
        })

    } else if(all(cd$data_type == "pb")){
        offset <- colMeans(sweep(bin_counts, 2, cd$ncells, "/"))
        cd$logitOffset <- log(offset/(1-offset))

        full_formula <- as.formula(paste("cbind(succes,failure)",
                                         paste(formula, collapse = " "),
                                         "+ offset(logitOffset)"))

        mod <- lapply(seq_len(nrow(bin_counts))[1:100], function(i){
            cd$succes <- bin_counts[i,]
            cd$failure <- cd$ncells - cd$succes
            mod <- glm(formula = full_formula,
                       family = "quasibinomial",
                       data = cd)
            return(mod)
        })
    }

    # test coefficient of interest
    sum_mod <- lapply(mod, summary)
    beta <- sapply(sum_mod, function(summod){coef(summod)[coef_test,1]})
    var_unscaled <- sapply(sum_mod, function(summod){summod$cov.unscaled[coef_test,2]})
    SE <- sapply(sum_mod, function(summod){coef(summod)[coef_test,2]})
    disp <- sapply(sum_mod, function(summod){summod$dispersion})
    pval <- sapply(sum_mod, function(summod){coef(summod)[coef_test,4]})
    res <- data.frame(gene = rownames(bin_counts)[1:100],
                      cluster_id = levels(cd$cluster_id),
                      logOdds = beta,
                      var_unscaled = var_unscaled,
                      dispersion = disp,
                      SE = SE,
                      Wald = (beta^2)/SE,
                      PValue = pval,
                      FDR = p.adjust(pval, method = "BH"),
                      comparsion = coef_test)
    return(res)
}


apply_qbGLM_offset_squeeze <- function(sce,
                                       coef_test = "group_idB",
                                       formula = ~ batch_cov + group_id,
                                       ...) {

    # Prepare for fit
    bin_counts <- as.matrix(assay(sce))
    cd <- colData(sce)

    if(all(cd$data_type == "sce")){
        offset <- colMeans(bin_counts)
        cd$logitOffset <- log(offset/(1-offset))

        formula <- as.formula(paste("y", paste(formula, collapse = " "),
                                    "+ offset(logitOffset)"))

        # fit models
        mod <- lapply(seq_len(nrow(bin_counts))[1:100], function(i){
            y <- bin_counts[i,]
            cd$y <- y
            mod <- glm(formula = formula, family = "quasibinomial", data = cd)
            return(mod)
        })

    } else if(all(cd$data_type == "pb")){
        offset <- colMeans(sweep(bin_counts, 2, cd$ncells, "/"))
        cd$logitOffset <- log(offset/(1-offset))

        full_formula <- as.formula(paste("cbind(succes,failure)",
                                         paste(formula, collapse = " "),
                                         "+ offset(logitOffset)"))

        mod <- lapply(seq_len(nrow(bin_counts))[1:100], function(i){
            cd$succes <- bin_counts[i,]
            cd$failure <- cd$ncells - cd$succes
            mod <- glm(formula = full_formula,
                       family = "quasibinomial",
                       data = cd)
            return(mod)
        })
    }

    # test coefficient of interest
    dfresid_mod <- sapply(mod, function(mod){mod$df.residual})
    sum_mod <- lapply(mod, summary)
    disp <- sapply(sum_mod, function(summod){summod$dispersion})
    beta <- sapply(sum_mod, function(summod){coef(summod)[coef_test,1]})
    var_unscaled <- sapply(sum_mod, function(summod){summod$cov.unscaled[coef_test,coef_test]})

    # squeeze variance
    hlp <- squeezeVar(var = disp,
                      df = dfresid_mod,
                      robust = FALSE)

    # Compute relevant statistics using squeezed variances
    df_post <- hlp$df.prior + dfresid_mod # posterior DF
    SE <- sqrt(var_unscaled * hlp$var.post) # SE based on squeezed dispersion
    t <- beta / SE # new t-stat
    pval <- pt(-abs(t), df_post) * 2 # new pval

    res <- data.frame(gene = rownames(bin_counts)[1:100],
                      cluster_id = levels(cd$cluster_id),
                      logOdds = beta,
                      var_unscaled = var_unscaled,
                      dispersion = disp,
                      dispersion_sq = hlp$var.post,
                      df_prior = hlp$df.prior,
                      df_resid = dfresid_mod,
                      SE = SE,
                      t = t,
                      Wald = (beta^2)/SE,
                      PValue = pval,
                      FDR = p.adjust(pval, method = "BH"),
                      comparsion = coef_test)
    return(res)
}
