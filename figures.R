here::i_am("figures.R")
out_dir <- here::here("figures")

# Setup ------------------------------------------------------------------------

library(DDCompanion)
library(iCOBRA)
library(SingleCellExperiment)
library(scuttle)

library(tidyverse)
library(ggthemes)
library(patchwork)

# Setup-------------------------------------------------------------------------

theme_set(theme_light(base_size = 14))
method_cols <- c(
    "edgeR_NB" = "#F28E2B",
    "edgeR_QP" = "red",
    "edgeR_NB_optim" = "gold1",
    "edgeR_QP_optim" = "hotpink3",
    "bGLM" = "grey32",
    "qbGLM" = "bisque3",
    "qbGLM_offset" = "dodgerblue",
    "qbGLM_offset_squeeze" = "dodgerblue4"
)
main_methods <- c("bGLM", "qbGLM", "qbGLM_offset", "qbGLM_offset_squeeze",
                  "edgeR_NB", "edgeR_QP", "edgeR_NB_optim", "edgeR_QP_optim")

# Obtain helper functions for generating the figures----------------------------
source("./figures/figures_helpers.R")
message("Loaded helper functions to genrate figures")

# Read and wrangle the data objects needed for the figures----------------------
message("Started reading and wrangling data for figures")
message("This will take a couple of minutes")
source("./figures/figures_prepare.R")
message("Finished reading and wrangling data for figures")

# Figure 1: Lupus T4_naive_10 mock results  ------------------------------------

fig_ncM_sub <- plot_mock(data = lupus_mock_table_ncM$ncM_10,
                          multisample = FALSE, stratifier_row = "method")
ggsave(file.path(out_dir, "fig1_lupus_ncM_mock.png"), fig_ncM_sub, width = 12, height=12*0.618)

# Figure S1-3: Lupus all cell types all sample size pseudobulk mock results-----
## We generate 3 plots, one for each cell type. Each plot plot has 8*4 panels,
## for 8 methods and 4 sample sizes.
fig_ncM <- plot_mock(data = as.data.frame(do.call(rbind, lupus_mock_table_ncM)),
                     multisample = TRUE, stratifier_col = "samplesize", stratifier_row = "method")
fig_T4naive <- plot_mock(data = as.data.frame(do.call(rbind, lupus_mock_table_T4naive)),
                         multisample = TRUE, stratifier_col = "samplesize", stratifier_row = "method")
fig_Bmem <- plot_mock(data = as.data.frame(do.call(rbind, lupus_mock_table_Bmem)),
                      multisample = TRUE, stratifier_col = "samplesize", stratifier_row = "method")

ggsave(file.path(out_dir, "figS1_lupus_ncM_mock.png"), fig_ncM, width = 12, height=28*0.618)
ggsave(file.path(out_dir, "figS2_lupus_T4naive_mock.png"), fig_T4naive, width = 12, height=28*0.618)
ggsave(file.path(out_dir, "figS3_lupus_Bmem_mock.png"), fig_Bmem, width = 12, height=28*0.618)

# Figure S4: Lupus all cell types, 22v22 single-cell data ----------------------
fig_sce_mock <- plot_mock(data = as.data.frame(do.call(rbind, lupus_mock_table_sce)),
                         multisample = TRUE, stratifier_col = "celltype", stratifier_row = "method")
ggsave(file.path(out_dir, "figS4_lupus_mock_sce.png"), fig_sce_mock, width = 12, height=28*0.618)

# Figure S5-7: Covid all celltype mock results------------------------------------
fig_covid_class_switched <- plot_mock(data = as.data.frame(do.call(rbind,
                                                                   list(covid_mock_tables$class_switched_memory_B_cell_Healthy,
                                                                        covid_mock_tables$class_switched_memory_B_cell_Mild,
                                                                        covid_mock_tables$class_switched_memory_B_cell_Moderate,
                                                                        covid_mock_tables$class_switched_memory_B_cell_Critical))),
                                      multisample = TRUE,
                                      stratifier_col = "disease",
                                      stratifier_row = "method")
ggsave(file.path(out_dir, "figS5_covid_class_switched_mock.png"), fig_covid_class_switched, width = 12, height=28*0.618)

fig_covid_immature <- plot_mock(data = as.data.frame(do.call(rbind,
                                                             list(covid_mock_tables$immature_B_cell_Healthy,
                                                                  covid_mock_tables$immature_B_cell_Mild,
                                                                  covid_mock_tables$immature_B_cell_Moderate,
                                                                  covid_mock_tables$immature_B_cell_Critical))),
                                multisample = TRUE,
                                stratifier_col = "disease",
                                stratifier_row = "method")
ggsave(file.path(out_dir, "figS6_covid_immature_mock.png"), fig_covid_immature, width = 12, height=28*0.618)

fig_covid_naive <- plot_mock(data = as.data.frame(do.call(rbind,
                                                          list(covid_mock_tables$naive_B_cell_Healthy,
                                                               covid_mock_tables$naive_B_cell_Mild,
                                                               covid_mock_tables$naive_B_cell_Moderate,
                                                               covid_mock_tables$naive_B_cell_Critical))),
                             multisample = TRUE,
                             stratifier_col = "disease",
                             stratifier_row = "method")
ggsave(file.path(out_dir, "figS7_covid_naive_mock.png"), fig_covid_naive, width = 12, height=28*0.618)

# Figure S9: Lupus sim all ct results ------------------------------------------
lupus_sim_cobra <- map2(lupus_sim_res, lupus_pb_objects,
    ~ map2(.x, .y, prepare_COBRAData))

## Supplementary figure
lupus_sim_perf <- map_depth(lupus_sim_cobra, 2,
    calculate_performance, binary_truth = "status")
lupus_sim_plot_objects <- map_depth(lupus_sim_perf, 2,
    prepare_data_for_plot, colorscheme = c(method_cols, truth = "grey35"))
lupus_sim_plot_data <- map(lupus_sim_plot_objects, combine_fdrtpr_tables) |>
    bind_rows(.id = "celltype")
lupus_sim_plot_data$replicate <- str_to_title(sub("_", " ", lupus_sim_plot_data$replicate))
lupus_sim_plot_data <- lupus_sim_plot_data[-which(lupus_sim_plot_data$method == "bGLM"),]

## plot
fig_lupus_sim_all <- lupus_sim_plot_data |>
    group_by(thr, method, celltype) |>
    summarize(across(c(FDR, TPR), mean), .groups = "keep") |>
    plot_fdrtpr_points() +
        base_theme(legend.position = "bottom") +
        facet_wrap(vars(factor(celltype, levels=c("ncM_10", "T4_naive_10", "B_mem_10",
                                                  "ncM_20", "T4_naive_20", "B_mem_20",
                                                  "ncM_30", "T4_naive_30", "B_mem_30",
                                                  "ncM_44", "T4_naive_44", "B_mem_44"))),
                   nrow = 4,
                   scales="fixed") +
    theme(strip.text = element_text(size=16))

## save
save_plot(file.path(out_dir, "figS9_lupus_sim_all.png"), fig_lupus_sim_all, width = 12, asp = 0.618*2)

# Figure 2: Lupus sim all ct 5v5 results ---------------------------------------
lupus_sim_plot_data <- map(lupus_sim_plot_objects, combine_fdrtpr_tables) |>
    bind_rows(.id = "celltype")
lupus_sim_plot_data$replicate <- str_to_title(sub("_", " ", lupus_sim_plot_data$replicate))
lupus_sim_plot_data <- lupus_sim_plot_data[lupus_sim_plot_data$celltype %in% c("ncM_10", "T4_naive_10", "B_mem_10"),]

## plot
fig_lupus_sim_5v5 <- lupus_sim_plot_data |>
    group_by(thr, method, celltype) |>
    summarize(across(c(FDR, TPR), mean), .groups = "keep") |>
    plot_fdrtpr_points() +
    base_theme(legend.position = "bottom") +
    facet_wrap(vars(factor(celltype, levels=c("ncM_10", "T4_naive_10", "B_mem_10"))),
               nrow = 1,
               scales="fixed") +
    theme(strip.text = element_text(size=16))

## save
save_plot(file.path(out_dir, "fig2_lupus_sim_5v5.png"), fig_lupus_sim_5v5, width = 12, asp = 0.618/1.33)

# Figure S8: Lupus simulation results all ct alternative -----------------------

use_thr <- 0.05
lupus_sim_plot_data_s9 <- lupus_sim_plot_data |>
    filter(thr == use_thr)

lupus_sim_plot_data_s9$method <- factor(lupus_sim_plot_data_s9$method,
                                levels = main_methods)

fig_lupus_sim_all_alt <- plot_fdr_control(lupus_sim_plot_data_s9, shape_by = "celltype") +
    scale_shape_manual(name = "Cell Type", values = 21:23)

save_plot(file.path(out_dir, "figS8_lupus_sim_all_alt.png"), fig_lupus_sim_all_alt, width = 12, asp = 0.618)

# Figure S10: Covid all celltypes sim results------------------------------------
covid_all_sim_res <- map(covid_all_sim_res,
                         ~ rename_methods(.x) |>
                             map(get_tables, depth = 1) |>
                             ## Bring replicates to top level
                             transpose() |>
                             map(combine_tables)
)
covid_cobra_data <- map2(covid_all_sim_res, covid_all_pb_objects,
                         ~ map2(.x, .y, prepare_COBRAData)
)
covid_cobra_perf <- map_depth(covid_cobra_data, 2,
                              calculate_performance, binary_truth = "status"
)
covid_cobra_objects <- map_depth(covid_cobra_perf, 2,
                                 prepare_data_for_plot, colorscheme = c(method_cols, truth = "grey35")
)

sim_covid_data <- map(covid_cobra_objects, combine_fdrtpr_tables) |>
    bind_rows(.id = "celltype")
sim_covid_data$replicate <- str_to_title(sub("_", " ", sim_covid_data$replicate))
sim_covid_data$celltype <- factor(sim_covid_data$celltype,
                                     levels = c("class_switched_memory_B_cell_Healthy",
                                                "class_switched_memory_B_cell_Mild",
                                                "class_switched_memory_B_cell_Moderate",
                                                "class_switched_memory_B_cell_Critical",
                                                "immature_B_cell_Healthy",
                                                "immature_B_cell_Mild",
                                                "immature_B_cell_Moderate",
                                                "immature_B_cell_Critical",
                                                "naive_B_cell_Healthy",
                                                "naive_B_cell_Mild",
                                                "naive_B_cell_Moderate",
                                                "naive_B_cell_Critical"))

fig_covid_sim_all <- sim_covid_data |>
    group_by(thr, method, celltype) |>
    summarize(across(c(FDR, TPR), mean), .groups = "keep") |>
    plot_fdrtpr_points() +
    base_theme(legend.position = "bottom") +
    facet_wrap(vars(celltype),
               ncol=4) +
    theme(strip.text = element_text(size=16))

save_plot(file.path(out_dir, "figS10_covid_sim_all.png"), fig_covid_sim_all, width = 22, asp = 0.618*1.33)
# TODO:  Maybe add the number of patients somewhere



