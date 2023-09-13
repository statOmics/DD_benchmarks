# Run a differential detection analysis on pseudobulk data using edgeR_NB_optim
DD_pb_edgeR_optim_NB <- function(object){

    # binarize the counts
    object_bin <- object
    assay(object_bin)[assay(object_bin)>=1] <- 1

    # aggregate
    pb_bin <- aggregateAcrossCells(object_bin, ids = object_bin$ind_cov)
    pb_bin$data_type <- "pb"

    # create formula
    formula <- ~ group_id

    # run DD analysis with edgeR_NB_optim
    res <- DDCompanion:::apply_edgeR_NB_optim(sce = pb_bin,
                                              coef_test = "group_idB",
                                              formula = formula)
    return(res)
}

# Run a differential expression analysis on pseudobulk data using edgeR
DE_pb_edgeR <- function(object){

    # aggregate
    pb <- aggregateAcrossCells(object, ids = object$ind_cov)

    # run DD analysis with edgeR_NB
    y <- DGEList(counts=assay(pb), samples=colData(pb))
    y <- calcNormFactors(y)
    design <- model.matrix(~ factor(batch_cov) + factor(mock_group), y$samples)
    y<-estimateDisp(y,design)
    fit <- glmQLFit(y,design,robust=T)

    res_DGE <- glmQLFTest(fit,coef = ncol(design))
    res_DGE$table$FDR <- p.adjust(res_DGE$table$PValue, method = "BH")

    return(res_DGE$table)
}

# Ad-hoc differential signal: strategy 1
# In the first ad-hoc strategy, for each gene in which signal will be introduced,
# we in 1 mock treatment group select the cells in the lowest 25% quantile that
# express the target gene, and set their expression value to zero. This strategy
# will generate DD signal, as the detection of the gene is reduced in 1 treatment
# group. While this strategy will also generate a DE signal, i.e. a shift in mean
# expression between the treatment groups, this signal will be relatively small,
# because only the lowest expression values are set to zero. We expect this signal
# to mimic real-life situations, in which genes that are truly lowly expressed in
# the cell are not detected by the transcriptomics technology.

DD_adhoc_1 <- function(object, genes_to_change, frac_cells){

    hlp_a <- assay(object)[genes_to_change[1:(length(genes_to_change)/2)],object$mock_group=="A"]
    hlp_b <- assay(object)[genes_to_change[(length(genes_to_change)/2+1):length(genes_to_change)],object$mock_group=="B"]

    vec_zero <- c(round(rowSums(hlp_a > 0) * frac_cells) ,
                  round(rowSums(hlp_b > 0) * frac_cells))

    for (i in 1:nrow(hlp_a)) {
        nonzero <- which(hlp_a[i,] > 0)
        n_zero <- round(length(nonzero)*frac_cells)
        count_vec <- sort(hlp_a[i,names(nonzero)])
        hlp_a[i,names(count_vec)[1:n_zero]] <- 0
    }

    for (i in 1:nrow(hlp_b)) {
        nonzero <- which(hlp_b[i,] > 0)
        n_zero <- round(length(nonzero)*frac_cells)
        count_vec <- sort(hlp_b[i,names(nonzero)])
        hlp_b[i,names(count_vec)[1:n_zero]] <- 0
    }

    assay(object)[rownames(hlp_a),object$mock_group=="A"] <- hlp_a
    assay(object)[rownames(hlp_b),object$mock_group=="B"] <- hlp_b

    rowData(object)$ncells_zero <- NA
    rowData(object[names(vec_zero),])$ncells_zero_2 <- vec_zero
    rowData(object[names(vec_zero),])$is_DE <- TRUE
    return(object)
}

## Ad-hoc differential signal: strategy 2
# The second ad-hoc strategy is similar to the one discussed above. However, the
# counts of the cells that are set to zero are uniformly distributed over the
# cells in the same treatment group that still express the gene. The effect on the
# DD signal is identical, but this strategy enforces that the within-group average
# expression of the target gene is not altered. However, note that this strategy
# will result in increased variability in the expression data, and may result in
# artificially generating bimodal distributions for target genes.

DD_adhoc_2 <- function(object, genes_to_change, frac_cells){

    hlp_a <- assay(object)[genes_to_change[1:(length(genes_to_change)/2)],object$mock_group=="A"]
    hlp_b <- assay(object)[genes_to_change[(length(genes_to_change)/2+1):length(genes_to_change)],object$mock_group=="B"]

    vec_zero <- c(round(rowSums(hlp_a > 0) * frac_cells) ,
                  round(rowSums(hlp_b > 0) * frac_cells))

    for (i in 1:nrow(hlp_a)) {
        nonzero <- names(which(hlp_a[i,] > 0))
        n_zero <- round(length(nonzero)*frac_cells)

        # change some to zero
        change_to_zero <- sample(nonzero, n_zero)
        former_counts <- sum(unname(hlp_a[i,change_to_zero])) # store sum of counts set to zero
        hlp_a[i,change_to_zero] <- 0

        # add back the removed counts to the other cells
        not_change_to_zero <- nonzero[!nonzero %in% change_to_zero]
        hlp_add <- table(sample(not_change_to_zero,
                                former_counts,
                                replace = TRUE)) # redistribute the formerly removed counts
        hlp_a[i,names(hlp_add)] <- as.matrix(hlp_a[i,names(hlp_add)] + hlp_add)
    }

    for (i in 1:nrow(hlp_b)) {
        nonzero <- names(which(hlp_b[i,] > 0))
        n_zero <- round(length(nonzero)*frac_cells)

        # change some to zero
        change_to_zero <- sample(nonzero, n_zero)
        former_counts <- sum(unname(hlp_b[i,change_to_zero])) # store sum of counts set to zero
        hlp_b[i,change_to_zero] <- 0

        # add back the removed counts to the other cells
        not_change_to_zero <- nonzero[!nonzero %in% change_to_zero]
        hlp_add <- table(sample(not_change_to_zero,
                                former_counts,
                                replace = TRUE)) # redistribute the formerly removed counts
        hlp_b[i,names(hlp_add)] <- as.matrix(hlp_b[i,names(hlp_add)] + hlp_add)
    }

    assay(object)[rownames(hlp_a),object$mock_group=="A"] <- hlp_a
    assay(object)[rownames(hlp_b),object$mock_group=="B"] <- hlp_b

    rowData(object)$ncells_zero_3 <- NA
    rowData(object[names(vec_zero),])$ncells_zero_3 <- vec_zero
    rowData(object[names(vec_zero),])$is_DE <- TRUE
    return(object)
}

## Ad-hoc differential signal: strategy 3
# In the third ad-hoc simulation strategy, we simply add a count of 1 to 25% of
# the cells in 1 treatment group that already expressed the gene. This way, we
# generate a DE signal without introducing DD.

DD_adhoc_3 <- function(object, genes_to_change, frac_cells){

    hlp_a <- assay(object)[genes_to_change[1:(length(genes_to_change)/2)],object$mock_group=="A"]
    hlp_b <- assay(object)[genes_to_change[(length(genes_to_change)/2+1):length(genes_to_change)],object$mock_group=="B"]

    vec_add <- c(round(rowSums(hlp_a > 0) * frac_cells) ,
                 round(rowSums(hlp_b > 0) * frac_cells))

    for (i in 1:nrow(hlp_a)) {
        nonzero <- which(hlp_a[i,] > 0)
        n_add <- round(length(nonzero)*frac_cells)
        target <- names(sample(nonzero, n_add))
        hlp_a[i,target] <- hlp_a[i,target] + 1
    }

    for (i in 1:nrow(hlp_b)) {
        nonzero <- which(hlp_b[i,] > 0)
        n_add <- round(length(nonzero)*frac_cells)
        target <- names(sample(nonzero, n_add))
        hlp_b[i,target] <- hlp_b[i,target] + 1
    }

    assay(object)[rownames(hlp_a),object$mock_group=="A"] <- hlp_a
    assay(object)[rownames(hlp_b),object$mock_group=="B"] <- hlp_b

    rowData(object)$ncells_add <- NA
    rowData(object[names(vec_add),])$ncells_add <- vec_add
    rowData(object[names(vec_add),])$is_DE <- TRUE
    return(object)
}

# Visualize results
visualize_stageR_fdrtpr <- function(DD_pval, DE_pval, truth){

    pvalues <- data.frame(DD_edgeR_NB_optim = DD_pval,
                          DE_edgeR_NB = DE_pval)
    #screening stage
    pvalues$stagewise <- apply(X = pvalues,
                               MARGIN = 1,
                               FUN = hmp.stat)
    rownames(pvalues) <- rownames(truth)
    cobra <- COBRAData(pval = pvalues, truth = truth)
    cobra <- calculate_adjp(cobra)
    cobraperf <- calculate_performance(cobra,
                                       binary_truth = "is_DE",
                                       splv = "none")
    cobraplot <- prepare_data_for_plot(cobraperf,
                                       colorscheme = "Dark2",
                                       facetted = TRUE)
    gg <- plot_fdrtprcurve(cobraplot) +
        theme(strip.background = element_blank(),
              strip.text.x = element_blank(),
              plot.title = element_text(size=10),
              axis.text.x = element_text(size = 10),
              axis.text.y = element_text(size = 10),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_text(size = 10))
    return(gg)
}









