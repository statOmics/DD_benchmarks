# helper functions to generate the figures

# Three functions making pretty plots
method_col_scale <- function(...) scale_color_manual(values = method_cols, ...)
method_fill_scale <- function(...) scale_fill_manual(values = method_cols, ...)
base_theme <- function(...) {
    theme(
        text = element_text(size = 16),
        # panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 11, face = "bold", color = "black"),
        strip.background = element_blank(),
        ...
    )
}

label_func_grid_row <- function(stratifier_row) {
    return(paste0(1:length(unique(stratifier_row)), ": ", stratifier_row))
}
label_func_grid_col <- function(stratifier_col) {
    return(paste0("Panel ", LETTERS[1:length(unique(stratifier_col))], ": ", stratifier_col))
}
label_func_wrap <- function(stratifier_row) {
    return(paste0("Panel ", LETTERS[1:length(unique(stratifier_row))], ": ", stratifier_row))
}

# TODO this is now a ghost function, but if I remove it downstream code seems to break...
rename_methods <- function(res_list) {
    res_list
}

# Function for saving plots
save_plot <- function(filename, plot, width = 9, asp = 0.618, ...) {
    height <- width * asp
    ggsave(filename, plot, width = width, height = height, ...)
}

# Function to concatenate fdrtpr tables across replicate datasets
combine_fdrtpr_tables <- function(cobra_list) {
    out <- map_dfr(cobra_list, fdrtpr, .id = "replicate")
    mutate(out,
           thr = as.numeric(sub("thr", "", thr)),
           method = as.factor(method)
    )
}

# Function to concatenate fdrtpr tables across replicate datasets
plot_fdrtpr_points <- function(x, thr = c(0.01, 0.05, 0.1)) {
    ggplot(x, aes(FDR, TPR, color = method, fill = method)) +
        geom_vline(xintercept = thr, linetype = "dashed", alpha = 0.5) +
        geom_path(linewidth = 0.8, alpha = 0.6, show.legend = FALSE) +
        geom_point(size = 4, stroke = 0.8, alpha = 0.8, shape = 21, col = "black") +
        labs(fill = NULL) +
        method_col_scale() +
        method_fill_scale()
}

# Concatenate FDR-TPR tables across replicate datasets
combine_fdrtpr_tables <- function(cobra_list) {
    out <- map_dfr(cobra_list, fdrtpr, .id = "replicate")
    mutate(out,
           thr = as.numeric(sub("thr", "", thr)),
           method = as.factor(method)
    )
}

# Make mock comparison
plot_mock <- function(data, multisample = c(TRUE, FALSE), stratifier_col, stratifier_row){
    if(multisample){
        data$stratifier_col <- data[[stratifier_col]]
        data$stratifier_row <- data[[stratifier_row]]
        ggplot(data = data,
               aes(PValue, col = stratifier_row)) +
            geom_density(aes(y = after_stat(scaled), lty = replicate),
                         adjust = 0.5, linewidth = 1, bounds = c(0,1),
                         key_glyph = "path"
            ) +
            facet_grid(rows = vars(stratifier_row),
                       cols = vars(stratifier_col),
                       labeller = labeller(stratifier_row = label_func_grid_row,
                                           stratifier_col = label_func_grid_col)) +
            method_col_scale() +
            method_fill_scale() +
            scale_x_continuous(breaks = seq(0, 1, 0.2)) +
            scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
            labs(x = "P-value", y = "Density") +
            guides(
                col = "none", fill = "none",
                lty = guide_legend(
                    title = NULL, nrow = 1,
                    override.aes = list(alpha = 1)
                )
            ) +
            base_theme(legend.position = "bottom")
    } else {
        data$stratifier_row <- data[[stratifier_row]]
        ggplot(data = data,
               aes(PValue, col = stratifier_row)) +
            geom_density(aes(y = after_stat(scaled), lty = replicate),
                         adjust = 0.5, linewidth = 1, bounds = c(0,1),
                         key_glyph = "path"
            ) +
            facet_wrap(stratifier_row,
                       nrow = 2,
                       labeller = labeller(.rows = label_func_wrap)) +
            method_col_scale() +
            method_fill_scale() +
            scale_x_continuous(breaks = seq(0, 1, 0.2)) +
            scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
            labs(x = "P-value", y = "Density") +
            guides(
                col = "none", fill = "none",
                lty = guide_legend(
                    title = NULL, nrow = 1,
                    override.aes = list(alpha = 1)
                )
            ) +
            base_theme(legend.position = "bottom")
    }

}


# Function that visualizes the working points when the FDR level is set at
# nominal levels of 1%, 5% and 10%, respectively, by drawing circles on the
# FDR-TPR curve
plot_fdrtpr_points <- function(x, thr = c(0.01, 0.05, 0.1)) {
    ggplot(x, aes(FDR, TPR, color = method, fill = method)) +
        geom_vline(xintercept = thr, linetype = "dashed", alpha = 0.5) +
        geom_path(linewidth = 0.8, alpha = 0.6, show.legend = FALSE) +
        geom_point(size = 4, stroke = 0.8, alpha = 0.8, shape = 21, col = "black") +
        labs(fill = NULL) +
        method_col_scale() +
        method_fill_scale()
}

# Function that visualizes the FDP and TPP of a performance analysis. In stead
# of generating a curve, it visualizes the raw FDP per dataset in a scatterplot,
# with the size of point representing the TPP
plot_fdr_control <- function(data, shape_by = NULL) {
    dat_text <- data.frame(
        label = paste("alpha == ", use_thr),
        method = levels(data$method)[1],
        celltype = levels(as.factor(data$celltype))[1],
        FDR = use_thr,
        TPP = 0
    )

    base_jitter <- function(...) {
        geom_jitter(
            width = 0.2, height = 0,
            alpha = 0.8, stroke = 0.8,
            ...
        )
    }
    if (!is.null(shape_by)) {
        shape_by <- sym(shape_by)
        add_jitter <- base_jitter(mapping = aes(shape = !!enquo(shape_by)))
    } else {
        add_jitter <- base_jitter(shape = 21)
    }
    data$TPP <- data$TPR
    ggplot(data, aes(celltype, FDR, fill = method, size = TPP)) +
        add_jitter +
        geom_hline(yintercept = use_thr, linetype = "dashed", alpha = 0.7) +
        facet_wrap(vars(factor(method,
                               levels = c("bGLM", "qbGLM",
                                          "qbGLM_offset", "qbGLM_offset_squeeze",
                                          "edgeR_NB", "edgeR_QP",
                                          "edgeR_NB_optim", "edgeR_QP_optim"))),
                   ncol = 4, scales = "free_x") +
        geom_text(
            data = dat_text, aes(label = label), size = 6,
            alpha = 0.5, hjust = 0.4, vjust = -0.5, parse = TRUE
        ) +
        labs(x = NULL, y = "FDP") +
        method_fill_scale(guide = "none") +
        theme_linedraw() +
        base_theme(
            panel.grid.major.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
        ) +
        theme(strip.background = element_rect(fill = "black"),
              strip.text = element_text(colour = "white"))
}

