#' Benchmark impact of similarity metric on clustering.
#'
#' @title scClustBench
#'
#' @examples
#' scClustBench(mat, method = "simlr", similarity = NULL, cores = 4, cores.ratio = 0)
#'
#' @param mat a (m x n) data matrix of gene expression measurements of individual cells with rows representing genes and columns representing cells. column names of the \emph{mat} must be cell types.
#' @param method Clustering method to be performed on the dataset between \emph{SIMLR} or \emph{kmeans} from \emph{amap} package. It is set to \emph{SIMLR} by default.
#' @param similarity A vector of similarity metrics to be used for clustering.
#' @param geneFilter A threshold to remove genes. The genes that are not expressed more than the threshold across all the cells in the dataset will be removed. Genes will not be removed if set to 0.
#' @param rep A number of subsampling of the matrix
#' @param subset_p Sampling percentage per cell types from \emph{mat}
#' @param cores Number of cores to be used for parallel processing
#' @param seed \emph{seed} for randomisation
#'
#' @return
#'
#'
#' @export scClustBench
#' @importFrom parallel stopCluster makeCluster detectCores clusterEvalQ
#' @importFrom parallel parLapply
#' @importFrom stats dnorm kmeans pbeta rnorm
#' @importFrom amap Kmeans
#' @importFrom methods is
#' @importFrom Hmisc rcorr
#' @importFrom caret createDataPartition
#' @importFrom mclust adjustedRandIndex
#' @importFrom clusteval cluster_similarity
#' @importFrom dendextend FM_index
#' @importFrom igraph compare
#' @import Matrix
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import doSNOW
#' @useDynLib scClust projsplx
#'
scClustBench <- function(mat, method = "simlr", similarity = NULL, geneFilter = 0.8, rep = 5, subset_p = 0.8, cores = 1, seed = 1, ...) {
  if (missing(mat)) {
    stop("No input for mat is given.")
  }

  if(Sys.info()[['sysname']]=="Windows"){
    os = "windows"
  } else {
    os = "other"
  }

  if (geneFilter != 0) {
    if (class(mat) == "data.frame") {
      mat <- as.matrix(mat)
    }
    if (typeof(mat) == "character") {
      mat <- toNumeric(mat)
    }

    # filter genes
    del <- rowSums(mat == 0) / ncol(mat) > geneFilter
    if (any(del == TRUE)) {
      mat <- mat[-which(del == TRUE),]
    }
  }

  mat <- filterMatrix(mat, gene_p = geneFilter)

  if (method == "simlr") {
    result <- simlrSubMatrix(mat, p = subset_p, similarity = similarity, rep = rep, cores = cores, seed = seed, os = os, ...)
  } else if (method == "kmeans") {
    if (similarity == "pearson") {
      similarity = "correlation"
    }
    result <- kmeanSubMatrix(mat, p = subset_p, similarity = similarity, rep = rep, cores = cores, seed = seed, os = os, ...)

  }
  return (result)
}
