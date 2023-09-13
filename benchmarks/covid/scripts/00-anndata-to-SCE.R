## Convert covid data to SCE and clean metadata
## ============================================

#renv::restore("../..")

message(getwd())
message(packageVersion("zellkonverter"))
message(packageVersion("CuratedAtlasQueryR"))
message(packageVersion("DDCompanion"))

## Directory setup
here_root <- "benchmarks/covid"
here::i_am(file.path(here_root, "scripts", "00-anndata-to-SCE.R"))

fs::dir_create(cache_dir <- here::here(here_root, "cache"))
fs::dir_create(out_dir <- here::here(here_root, "data"))

out_file <- file.path(out_dir, "covid-SCE-cleaned.rds")

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(CuratedAtlasQueryR)
    library(dplyr)
})

# Get SingleCellExperiment for Covid data -------------------------------------
metadata <- get_metadata()
metadata_target <- metadata |>
    dplyr::filter(
        name == "Single-cell multi-omics analysis of the immune response in COVID-19",
        cell_type %in% c("B cell",
                         "naive B cell",
                         "immature B cell",
                         "unswitched memory B cell",
                         "class switched memory B cell"),
        disease %in% c("COVID-19",
                       "normal")
    )
sce <- metadata_target |>
    get_single_cell_experiment(cache_directory = cache_dir)
colnames(colData(sce))[2] <- "cell_type_curated"

# Get author categories metadata ----------------------------------------
unharmonised_meta <- get_unharmonised_metadata(metadata_target,
                                               cache_directory = cache_dir)
cd_authors <- as.data.frame(unharmonised_meta$unharmonised)

relevant_author_vars <- c("pct_counts_mt", "initial_clustering", "Resample",
                          "Swab_result", "Status", "Smoker",
                          "Status_on_day_collection", "Collection_Day",
                          "Status_on_day_collection_summary", "Days_from_onset",
                          "Site", "Worst_Clinical_Status", "Outcome",
                          "donor_id", "author_cell_type")
cd_authors <- cd_authors[match(sce$original_cell_id, cd_authors$cell_),
                         relevant_author_vars]
colData(sce) <- cbind(colData(sce),
                      cd_authors)

# Filtering on author coldata -----------------------------------------------
sce <- sce[,sce$Status %in% c("Covid", "Healthy")]
colData(sce) <- droplevels(colData(sce))

# Export data -------------------------------------------------------------

## Save data as HDF5-backed SCE. This will save the object with pointers to the
## HDF5 file from which the object is derived.
## See also https://github.com/bioconductor/tenxbraindata#saving-computations

# move away from H5AD backend and save as regular rds
counts(sce) <- as(counts(sce), "sparseMatrix")
saveRDS(sce, out_file)

## "Reading this file back in with readRDS should work as long as the HDF5 file
## remains in its original location"
