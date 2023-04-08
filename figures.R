here::i_am("figures.R")
out_dir <- here::here("figures")
fs::dir_create(out_dir)

# Setup -------------------------------------------------------------------

library(DDCompanion)
library(iCOBRA)
library(SingleCellExperiment)

library(tidyverse)
library(ggthemes)
library(patchwork)
theme_set(theme_light(base_size = 14))

## Set color palette for methods
method_cols <- c(
    "edgeR_QP" = "#4E79A7",
    "edgeR_NB" = "#F28E2B",
    "bGLM" = "#E15759",
    "qbGLM" = "#76B7B2",
    "qbGLM_offset" = "#59A14F",
    "qbGLM_offset_squeeze" = "black"
)
method_col_scale <- function(...) scale_color_manual(values = method_cols, ...)
method_fill_scale <- function(...) scale_fill_manual(values = method_cols, ...)

main_methods <- c("edgeR_QP", "edgeR_NB", "bGLM", "qbGLM", "qbGLM_offset", "qbGLM_offset_squeeze")
pie(rep(1, length(main_methods)), labels = main_methods, col = method_cols)

base_theme <- function(...) {
    theme(
        text = element_text(size = 16),
        # panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 10, face = "bold", color = "black"),
        strip.background = element_blank(),
        ...
    )
}

## TODO this is now a ghost function, but if I remove it downstream code seems to break...
rename_methods <- function(res_list) {
    res_list
}

save_plot <- function(filename, plot, width = 9, asp = 0.618, ...) {
   height <- width * asp
   ggsave(filename, plot, width = width, height = height, ...)
}

# Load results ------------------------------------------------------------

## Lupus B_mem mock results
lupus_B_mem_mock_res <- get_mock_res_files(
    "lupus", main_methods, celltype = "B_mem", datatype = "pb"
) |>
    map(readRDS)

## Simulation results; use prop_DE = 0.05 for main figures
prop_DE <- 0.05

## Lupus B_mem sim results
lupus_B_mem_sim_res <- get_sim_res_files("lupus",
    main_methods, prop_DE = prop_DE, celltype = "B_mem", datatype = "pb"
) |>
    map(readRDS)

## All Lupus celltypes together
use_ct <- c("B_mem", "ncM") #  "T4_naive"
lupus_all_sim_res <- map(use_ct,
    ~ get_sim_res_files(
        dataset = "lupus", methods = main_methods,
        prop_DE = prop_DE, celltype = .x, datatype = "pb"
    ) |> map(readRDS)
) |>
    set_names(use_ct)

## All Lupus celltypes together
lupus_all_sce_objects <- map(use_ct,
    ~ get_SCE_files(
        dataset = "lupus", which = "sim_replicates",
        celltype = .x, prop_DE = prop_DE,
    )
) |>
    map(readRDS) |>
    set_names(use_ct)

# Figure 1: Mock results -------------------------------------------------

fig1_res <- lupus_B_mem_mock_res[main_methods] # |> rename_methods()
fig1_table <- map(fig1_res, get_aggregated_rep_tables, depth = 1) |>
    combine_tables()
fig1_table$replicate <- str_to_title(sub("_", " ", fig1_table$replicate))

fig1 <- ggplot(fig1_table, aes(PValue, col = method, fill = method)) +
    geom_density(aes(y = after_stat(scaled), lty = replicate),
        adjust = 0.5, size = 1, alpha = 0.1,
        key_glyph = "path"
    ) +
    # facet_grid(cols = vars(method), rows = vars(dataset)) +
    facet_grid(cols = vars(method)) +
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

ggsave(file.path(out_dir, "figure1.pdf"), fig1, width = 12, height=12*0.618/2)

# Figure 2: Simulation results -----------------------------------

lupus_all_sim_res <- map(lupus_all_sim_res,
                         ~ rename_methods(.x) |>
                             map(get_tables, depth = 1) |>
                             ## Bring replicates to top level
                             transpose() |>
                             map(combine_tables)
)

lupus_cobra_data <- map2(lupus_all_sim_res, lupus_all_sce_objects,
    ~ map2(.x, .y, prepare_COBRAData)
)
lupus_cobra_perf <- map_depth(lupus_cobra_data, 2,
    calculate_performance, binary_truth = "status"
)
lupus_cobra_objects <- map_depth(lupus_cobra_perf, 2,
    prepare_data_for_plot, colorscheme = c(method_cols, truth = "grey35")
)

combine_fdrtpr_tables <- function(cobra_list) {
    out <- map_dfr(cobra_list, fdrtpr, .id = "replicate")
    mutate(out,
        thr = as.numeric(sub("thr", "", thr)),
        method = as.factor(method)
    )
}

plot_fdrtpr_points <- function(x, thr = c(0.01, 0.05, 0.1)) {
    ggplot(x, aes(FDR, TPR, color = method, fill = method)) +
        geom_vline(xintercept = thr, linetype = "dashed", alpha = 0.5) +
        geom_path(size = 0.8, alpha = 0.6, show.legend = FALSE) +
        geom_point(size = 4, stroke = 0.8, alpha = 0.8, shape = 21, col = "black") +
        labs(fill = NULL) +
        method_col_scale() +
        method_fill_scale()
}

plot_fdr_control <- function(data, shape_by = NULL) {
    dat_text <- data.frame(
        label = paste("alpha == ", use_thr),
        method = levels(data$method)[1],
        FDR = use_thr, TPR = 0
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

    ggplot(data, aes(method, FDR, fill = method, size = TPR)) +
        add_jitter +
        geom_hline(yintercept = use_thr, linetype = "dashed", alpha = 0.7) +
        facet_wrap(vars(method), ncol = nlevels(data$method), scales = "free_x") +
        geom_text(
            data = dat_text, aes(label = label), size = 6,
            alpha = 0.5, hjust = 1, vjust = -0.5, parse = TRUE
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


fig2_lupus_data <- map(lupus_cobra_objects, combine_fdrtpr_tables) |>
    bind_rows(.id = "celltype")
fig2_lupus_data$replicate <- str_to_title(sub("_", " ", fig2_lupus_data$replicate))

fig2a <- fig2_lupus_data |>
    group_by(thr, method, celltype) |>
    summarize(across(c(FDR, TPR), mean), .groups = "keep") |>
    plot_fdrtpr_points() +
        base_theme(legend.position = "bottom") +
        facet_wrap(vars(celltype))

use_thr <- 0.05
fdr_lupus_data <- fig2_lupus_data |>
    filter(thr == use_thr)

fig2b <- plot_fdr_control(fdr_lupus_data, shape_by = "celltype") +
    scale_shape_manual(name = "Cell Type", values = 21:23)

fig2 <- wrap_plots(A = fig2a, B = fig2b, ncol = 2, guides = "collect") +
    plot_annotation(theme = theme(legend.position = "bottom"),
                    tag_levels = list(c("A", "B")))

save_plot(file.path(out_dir, "figure2.pdf"), fig2, width = 16, asp = 0.618/2)

# Figure 3: Influence # of subjects ----------------------------

n_patients <- c(10, 20, 30)
fig3_res_list <- map(n_patients,
    ~ get_sim_res_files(
        dataset = "lupus-n_patients",
        methods = main_methods,
        prop_DE = prop_DE,
        celltype = "ncM",
        datatype = "pb",
        n_patients = .x
    ) |> rename_methods()
) |>
    set_names(paste0("n_patients_", n_patients)) |>
    map_depth(2, readRDS)

fig3_res_tables <- map_depth(fig3_res_list, 2, get_aggregated_rep_tables, depth = 1) |>
    map(combine_tables, .id = "method")

fig3_sce_objects <- map(n_patients, ~ get_SCE_files(
    dataset = "lupus-n_patients", which = "sim_replicates",
    celltype = "ncM", n_patients = .x, prop_DE = prop_DE
)) |>
    set_names(paste0("n_patients_", n_patients)) |>
    map(readRDS)

fig3_cobra_data <- map2(fig3_res_tables, fig3_sce_objects, function(res_table, sce_list) {
    ## Split up results per replicate
    res_per_replicate <- split(res_table, res_table$replicate)
    map2(res_per_replicate, sce_list, prepare_COBRAData, replace_missing = TRUE)
})

fig3_cobra_perf <- map_depth(fig3_cobra_data, 2, calculate_performance, binary_truth = "status")
fig3_cobra_objects <- map_depth(fig3_cobra_perf, 2, prepare_data_for_plot)

## Add the original results (all patients) as well
fig3_cobra_objects[["n_patients_44"]] <- lupus_cobra_objects$ncM


## ----fig4, fig.asp=0.618/2-------------------------------------------------------
fig3_fdr_tpr_points <- map(fig3_cobra_objects, combine_fdrtpr_tables) |>
    bind_rows(.id = "n_patients") |>
    mutate(n_patients = as.numeric(sub("n_patients_", "", n_patients)))

fig3_fdr_tpr_points_averaged <- fig3_fdr_tpr_points |>
    group_by(thr, method, n_patients) |>
    summarize(across(c(FDR, TPR), mean), .groups = "keep")

fig3 <- plot_fdrtpr_points(fig3_fdr_tpr_points_averaged) +
    facet_wrap(vars(n_patients), nrow = 1, labeller = function(x) label_both(x, sep = " = ")) +
    scale_x_continuous(limits = c(0, 0.42)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
    base_theme() +
    theme(strip.background = element_rect(fill = "black"),
          strip.text = element_text(colour = "white"))

save_plot(file.path(out_dir, "figure3.pdf"), fig3, width = 12, asp = 0.618/2)
