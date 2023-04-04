#' Mock up a multi-subject SingleCellExperiment
#'
#' Calls [`scuttle::mockSCE()`] and adds subject and group labels.
#'
#' @param ncells The desired total number of cells.
#' @param ngenes The desired total number of genes.
#' @param nsubjects The desired total number of subjects.
#' @param ngroups The desired total number of **subject-level** groups.
#' @param nbatches The desired total number of **subject-level** batches.
#' @param ... Further arguments passed on to [`scuttle::mockSCE()`].
#'
#' @return
#' A \linkS4class{SingleCellExperiment} object. Subject, group and batch labels
#' are stored in the `"subject_id"`, `"mock_group"` and `"batch_id"` columns of
#' the `colData`, respectively. Subjects will be equally divided among the
#' groups.
#' @export
mockMultiSubjectSCE <- function(ncells = 200, ngenes = 2000,
                    nsubjects = 8, ngroups = 2, nbatches = 1,
                    ...) {
    out <- scuttle::mockSCE(ncells = ncells, ngenes = ngenes, ...)
    ## Remove 'Treatment' column to avoid confusion
    out$Treatment <- NULL

    out$subject_id <- gl(nsubjects, ncells / nsubjects,
        labels = sprintf("subject_%i", seq_len(nsubjects))
    )

    ## Assign batches
    batches <- gl(nbatches, nsubjects / nbatches,
        labels = sprintf("batch_%i", seq_len(nbatches))
    )
    out$batch_id <- batches[out$subject_id]

    ## Assign mock groups, dividing subjects in an unpaired design
    createMockGroups(out,
        labs = LETTERS[seq_len(ngroups)],
        subject_id = "subject_id", batch_id = "batch_id"
    )
}
