#' performs clustering on single-cell RNA-seq data.
#'
#' @title scClust
#'
#' @examples
#' scClust(mat = BuettnerFlorian$in_X, c = BuettnerFlorian$n_clust, method = "simlr", similarity = "pearson")
#'
#' @param mat a (m x n) data matrix of gene expression measurements of individual cells with rows representing genes and columns representing cells.
#' @param nCs number of clusters to be estimated
#' @param method Clustering method to be performed on the dataset between \emph{SIMLR} or \emph{kmeans} from \emph{amap} package. It is set to \emph{SIMLR} by default.
#' @param similarity A similarity metric to be used for clustering. This must be one of "euclidean", "pearson" or "spearman". If \emph{kmeans} is chosen for method, other metrics available to kmeans in amap can be used.
#' @param geneFilter A threshold to remove genes. The genes that are not expressed more than the threshold across all the cells in the dataset will be removed. Genes will not be removed if set to 0.
#'
#' @return For SIMLR, list of 8 elements describing the clusters obtained by SIMLR, of which y are the resulting clusters:
#'		y = results of k-means clusterings,
#'  	S = similarities computed by SIMLR,
#'  	F = results from network diffiusion,
#'  	ydata = data referring the the results by k-means,
#'  	alphaK = clustering coefficients,
#'  	execution.time = execution time of the present run,
#'  	converge = iterative convergence values by T-SNE,
#'  	LF = parameters of the clustering
#'
#' @return For Kmeans, a list with components:
#'    cluster  = A vector of integers indicating the cluster to which each point is allocated
#'    centers  = A matrix of cluster centres
#'    withinss = The within-cluster sum of squares distance for each cluster
#'    size     = The number of points in each cluster
#'
#' @export scClust
#' @importFrom parallel stopCluster makeCluster detectCores clusterEvalQ
#' @importFrom parallel parLapply
#' @importFrom stats dnorm kmeans pbeta rnorm
#' @importFrom amap Kmeans
#' @importFrom methods is
#' @importFrom Hmisc rcorr
#' @import Matrix
#' @useDynLib scClust projsplx
#'
scClust <- function(mat, nCs, method = "simlr", similarity = "pearson", geneFilter = 0.8, seed = 1, ...) {
  if (missing(mat)) {
    stop("No input for mat is given.")
  }
  if (missing(nCs)) {
    stop("No input for k is given")
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

  if (method == "simlr") {
    set.seed(seed)
    result = SIMLR(mat, nCs, distMethod = similarity, ...)
  } else if (method == "kmeans") {
    if (similarity == "pearson") {
      similarity = "correlation"
    }
    set.seed(seed)
    result <- amap::Kmeans(t(mat), nCs, method = similarity, ...)
  }

  return(result)


}
