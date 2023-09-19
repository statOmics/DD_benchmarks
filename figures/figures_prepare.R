# Script to prepare read and wrangle the objects needed for the paper figures

# Setup -------------------------------------------------------------------

# here::i_am("figures_prepare.R")
# out_dir <- here::here("figures")
# fs::dir_create(out_dir)
#
# main_methods <- c("bGLM", "qbGLM", "qbGLM_offset", "qbGLM_offset_squeeze",
#                   "edgeR_NB", "edgeR_QP", "edgeR_NB_optim", "edgeR_QP_optim")

# Load Lupus results -----------------------------------------------------------

## Lupus all 44 patients: mock results for all cell types, pseudobulk
use_ct_lupus <- c("T4_naive", "B_mem", "ncM")
lupus_44_pb_mock_res <- map(use_ct_lupus,
                          ~ get_mock_res_files(
                              dataset = "lupus", methods = main_methods,
                              celltype = .x, datatype = "pb"
                          ) |> map(readRDS)
) |>
    set_names(paste0(use_ct_lupus, "_44"))

## Lupus all 44 patients: mock results for all cell types, single-cell
lupus_44_sce_mock_res <- map(use_ct_lupus,
                            ~ get_mock_res_files(
                                dataset = "lupus", methods = main_methods,
                                celltype = .x, datatype = "sce"
                            ) |> map(readRDS)
) |>
    set_names(paste0(use_ct_lupus, "_44"))

## Lupus all 44 patients: simulation results for all cell types, pseudobulk
lupus_44_sim_res <- map(use_ct_lupus,
                         ~ get_sim_res_files(
                             dataset = "lupus", methods = main_methods,
                             prop_DE = 0.05, celltype = .x, datatype = "pb"
                         ) |> map(readRDS)
) |>
    set_names(paste0(use_ct_lupus, "_44"))

## Lupus all 44 patients: pseudobulk datasets for all cell types
lupus_44_sce_objects <- map(use_ct_lupus,
                             ~ get_SCE_files(
                                 dataset = "lupus",
                                 which = "sim_replicates",
                                 celltype = .x,
                                 prop_DE = 0.05,
                             )
) |>
    map(readRDS) |>
    set_names(paste0(use_ct_lupus, "_44"))
lupus_44_pb_objects <- lapply(lupus_44_sce_objects, function(lupus_ct){
    lapply(lupus_ct, function(lupus_ct_rep){
        pb <- aggregateAcrossCells(lupus_ct_rep,
                                   ids = lupus_ct_rep$ind_cov)
        colnames(pb) <- paste0("identifier_", 1:ncol(pb))
        pb$data_type <- "pb"
        return(pb)
    })
})
names(lupus_44_pb_objects) <- names(lupus_44_sce_objects)
rm(lupus_44_sce_objects)

# Load Lupus n_patients results-------------------------------------------------

n_patients <- c(10, 20, 30, 10, 20, 30, 10, 20, 30)
use_ct_lupus_npatients <- c("T4_naive", "T4_naive", "T4_naive",
                            "ncM", "ncM", "ncM",
                            "B_mem", "B_mem", "B_mem")

## Lupus 10-30 patients: mock results for all cell types
lupus_npatients_mock_res <- map2(use_ct_lupus_npatients, n_patients, ~ get_mock_res_files(
    dataset = "lupus-n_patients",
    methods = main_methods,
    datatype = "pb",
    celltype = .x,
    n_patients = .y) |>
        map(readRDS))|>
    set_names(paste0(use_ct_lupus_npatients, "_", n_patients))

## Lupus 10-30 patients: simulation results for all cell types
lupus_npatients_sim_res <- map2(use_ct_lupus_npatients, n_patients, ~ get_sim_res_files(
    dataset = "lupus-n_patients",
    methods = main_methods,
    datatype = "pb",
    prop_DE = 0.05,
    celltype = .x,
    n_patients = .y) |>
        map(readRDS))|>
    set_names(paste0(use_ct_lupus_npatients, "_", n_patients))

## Lupus 10-30 patients: pseudobulk datasets for all cell types
lupus_npatients_all_pb_objects <- map2(use_ct_lupus_npatients, n_patients, ~ get_SCE_files(
    dataset = "lupus-n_patients",
    which = "sim_replicates",
    celltype = .x,
    n_patients = .y,
    prop_DE = 0.05
)) |>
    set_names(paste0(use_ct_lupus_npatients, "_", n_patients)) |>
    map(readRDS)

# Combine Lupus and lupus-n_patients mock results per cell type ----------------

## T4naive
lupus_mock_res_T4naive <- list(T4_naive_10 = lupus_npatients_mock_res$T4_naive_10,
                               T4_naive_20 = lupus_npatients_mock_res$T4_naive_20,
                               T4_naive_30 = lupus_npatients_mock_res$T4_naive_30,
                               T4_naive_44 = lupus_44_pb_mock_res$T4_naive_44)
lupus_mock_table_T4naive <- lapply(seq_along(lupus_mock_res_T4naive), function(i,main_methods,samplesize){
    fig_i_res <- lupus_mock_res_T4naive[[i]][main_methods] # |> rename_methods()
    fig_i_table <- map(fig_i_res, get_aggregated_rep_tables, depth = 1) |>
        combine_tables()
    fig_i_table$replicate <- str_to_title(sub("_", " ", fig_i_table$replicate))
    fig_i_table$method <- factor(fig_i_table$method,
                                 levels = main_methods)
    fig_i_table$samplesize <- as.factor(samplesize[i])
    return(fig_i_table)
},main_methods=main_methods,samplesize=c("n = 10","n = 20","n = 30","n = 44"))
names(lupus_mock_table_T4naive) <- c("T4naive_10","T4naive_20","T4naive_30","T4naive_44")

## ncM
lupus_mock_res_ncM <- list(ncM_10 = lupus_npatients_mock_res$ncM_10,
                           ncM_20 = lupus_npatients_mock_res$ncM_20,
                           ncM_30 = lupus_npatients_mock_res$ncM_30,
                           ncM_44 = lupus_44_pb_mock_res$ncM_44)
lupus_mock_table_ncM <- lapply(seq_along(lupus_mock_res_ncM), function(i,main_methods,samplesize){
    fig_i_res <- lupus_mock_res_ncM[[i]][main_methods]
    fig_i_table <- map(fig_i_res, get_aggregated_rep_tables, depth = 1) |>
        combine_tables()
    fig_i_table$replicate <- str_to_title(sub("_", " ", fig_i_table$replicate))
    fig_i_table$method <- factor(fig_i_table$method,
                                 levels = main_methods)
    fig_i_table$samplesize <- as.factor(samplesize[i])
    return(fig_i_table)
},main_methods=main_methods,samplesize=c("n = 10","n = 20","n = 30","n = 44"))
names(lupus_mock_table_ncM) <- c("ncM_10","ncM_20","ncM_30","ncM_44")

## B_mem
lupus_mock_res_Bmem <- list(B_mem_10 = lupus_npatients_mock_res$B_mem_10,
                            B_mem_20 = lupus_npatients_mock_res$B_mem_20,
                            B_mem_30 = lupus_npatients_mock_res$B_mem_30,
                            B_mem_44 = lupus_44_pb_mock_res$B_mem_44)
lupus_mock_table_Bmem <- lapply(seq_along(lupus_mock_res_Bmem), function(i,main_methods,samplesize){
    fig_i_res <- lupus_mock_res_Bmem[[i]][main_methods]
    fig_i_table <- map(fig_i_res, get_aggregated_rep_tables, depth = 1) |>
        combine_tables()
    fig_i_table$replicate <- str_to_title(sub("_", " ", fig_i_table$replicate))
    fig_i_table$method <- factor(fig_i_table$method,
                                 levels = main_methods)
    fig_i_table$samplesize <- as.factor(samplesize[i])
    return(fig_i_table)
},main_methods=main_methods,samplesize=c("n = 10","n = 20","n = 30","n = 44"))
names(lupus_mock_table_Bmem) <- c("Bmem_10","Bmem_20","Bmem_30","Bmem_44")
rm(lupus_npatients_mock_res, lupus_44_pb_mock_res, lupus_mock_res_T4naive,
   lupus_mock_res_ncM,lupus_mock_res_Bmem)

lupus_sim_res <- list(T4_naive_10 = lupus_npatients_sim_res$T4_naive_10,
                      T4_naive_20 = lupus_npatients_sim_res$T4_naive_20,
                      T4_naive_30 = lupus_npatients_sim_res$T4_naive_30,
                      T4_naive_44 = lupus_44_sim_res$T4_naive_44,
                      ncM_10 = lupus_npatients_sim_res$ncM_10,
                      ncM_20 = lupus_npatients_sim_res$ncM_20,
                      ncM_30 = lupus_npatients_sim_res$ncM_30,
                      ncM_44 = lupus_44_sim_res$ncM_44,
                      B_mem_10 = lupus_npatients_sim_res$B_mem_10,
                      B_mem_20 = lupus_npatients_sim_res$B_mem_20,
                      B_mem_30 = lupus_npatients_sim_res$B_mem_30,
                      B_mem_44 = lupus_44_sim_res$B_mem_44)
lupus_sim_res <- map(lupus_sim_res,
                     ~ rename_methods(.x) |>
                         map(get_tables, depth = 1) |>
                         ## Bring replicates to top level
                         transpose() |>
                         map(combine_tables))
rm(lupus_npatients_sim_res, lupus_44_sim_res)

lupus_pb_objects <- list(T4_naive_10 = lupus_npatients_all_pb_objects$T4_naive_10,
                         T4_naive_20 = lupus_npatients_all_pb_objects$T4_naive_20,
                         T4_naive_30 = lupus_npatients_all_pb_objects$T4_naive_30,
                         T4_naive_44 = lupus_44_pb_objects$T4_naive_44,
                         ncM_10 = lupus_npatients_all_pb_objects$ncM_10,
                         ncM_20 = lupus_npatients_all_pb_objects$ncM_20,
                         ncM_30 = lupus_npatients_all_pb_objects$ncM_30,
                         ncM_44 = lupus_44_pb_objects$ncM_44,
                         B_mem_10 = lupus_npatients_all_pb_objects$B_mem_10,
                         B_mem_20 = lupus_npatients_all_pb_objects$B_mem_20,
                         B_mem_30 = lupus_npatients_all_pb_objects$B_mem_30,
                         B_mem_44 = lupus_44_pb_objects$B_mem_44)
rm(lupus_npatients_all_pb_objects, lupus_44_pb_objects)

# Mock lupus results on the single-cell level
lupus_44_sce_mock_res <- lupus_44_sce_mock_res[c(3,1,2)]
lupus_mock_table_sce <- lapply(seq_along(lupus_44_sce_mock_res), function(i,main_methods,celltype){
    fig_i_res <- lupus_44_sce_mock_res[[i]][main_methods]
    fig_i_table <- map(fig_i_res, get_aggregated_rep_tables, depth = 1) |>
        combine_tables()
    fig_i_table$replicate <- str_to_title(sub("_", " ", fig_i_table$replicate))
    fig_i_table$method <- factor(fig_i_table$method,
                                 levels = main_methods)
    fig_i_table$celltype <- as.factor(celltype[i])
    return(fig_i_table)
},main_methods=main_methods,celltype=c("ncM","T4_naive","B_mem"))
names(lupus_mock_table_sce) <- names(lupus_44_sce_mock_res)
rm(lupus_44_sce_mock_res)

# Load Covid results -----------------------------------------------------------

## All Covid mock results
use_ct_covid <- c("class_switched_memory_B_cell_Healthy",
                  "immature_B_cell_Healthy",
                  "naive_B_cell_Healthy",
                  "class_switched_memory_B_cell_Mild",
                  "immature_B_cell_Mild",
                  "naive_B_cell_Mild",
                  "class_switched_memory_B_cell_Moderate",
                  "immature_B_cell_Moderate",
                  "naive_B_cell_Moderate",
                  "class_switched_memory_B_cell_Critical",
                  "immature_B_cell_Critical",
                  "naive_B_cell_Critical")
covid_all_mock_res <- map(use_ct_covid,
                          ~ get_mock_res_files(
                              dataset = "covid", methods = main_methods,
                              celltype = .x, datatype = "pb"
                          ) |> map(readRDS)
) |>
    set_names(use_ct_covid)

covid_mock_tables <- lapply(seq_along(covid_all_mock_res), function(i,main_methods,celltype,disease){
    fig_i_res <- covid_all_mock_res[[i]][main_methods]
    fig_i_table <- map(fig_i_res, get_aggregated_rep_tables, depth = 1) |>
        combine_tables()
    fig_i_table$replicate <- str_to_title(sub("_", " ", fig_i_table$replicate))
    fig_i_table$method <- factor(fig_i_table$method,
                                 levels = main_methods)
    fig_i_table$celltype <- as.factor(celltype[i])
    fig_i_table$disease <- as.factor(disease[i])
    return(fig_i_table)
},main_methods=main_methods,
  celltype=rep(c("class_switched_memory", "immature", "naive"),times=4),
  disease = rep(c("Healthy","Mild","Moderate","Critical"),each=3)
)
names(covid_mock_tables) <- names(covid_all_mock_res)

## All Covid simulation results
covid_all_sim_res <- map(use_ct_covid,
                         ~ get_sim_res_files(
                             dataset = "covid", methods = main_methods,
                             prop_DE = 0.05, celltype = .x, datatype = "pb"
                         ) |> map(readRDS)
) |>
    set_names(use_ct_covid)

## All Covid pseudobulk datasets
covid_all_pb_objects <- map(use_ct_covid,
                             ~ get_SCE_files(
                                 dataset = "covid", which = "sim_replicates",
                                 celltype = .x, prop_DE = 0.05,
                             )
) |>
    map(readRDS) |>
    set_names(use_ct_covid)

gc()
