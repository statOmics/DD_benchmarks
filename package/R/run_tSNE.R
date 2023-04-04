#' @export
#' @importFrom scuttle logNormCounts
#' @importFrom scran modelGeneVarByPoisson getTopHVGs denoisePCA
#' @importFrom scater runTSNE
compute_tSNE <- function(x) {
    x <- logNormCounts(x)
    dec <- modelGeneVarByPoisson(x)
    top <- getTopHVGs(dec, prop = 0.1)
    x <- denoisePCA(x, technical = dec, subset.row = top)
    runTSNE(x, dimred = "PCA")
}
