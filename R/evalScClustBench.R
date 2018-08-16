#' Evaluate kmeans cluster result from \emph{scClustBench}
#'
#' @title evaluateScClustBench
#'
#' @examples
#' evaluateScClustBench(sub.simlr, method = "simlr")
#' evaluateScClustBench(sub.kmeans, method = "kmeans")
#'
#' @param sub.kmeans An object returned from \emph{scClustBench}
#' @param method A method used in \emph{scClustBench}
#'
#' @description Evaluate cluster results from \emph{scClustBench} with evaluation metrics NMI (Normalised Mutual Information), ARI (Adjusted Rand Index), FM (Fowlkes–Mallows) index and Jaccard Index
#'
#' @return A data frame with cluster evaluation for each similarity metric and for each repeats.
#'
#' @export evaluateScClustBench
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
"evalScClustBench" <- function(sub.kmeans, method = "simlr") {
  # for each subset matrix kmeans
  result.dat <- data.frame()
  for (i in 1:length(sub.kmeans)) {
    ari <- getARI(sub.kmeans[[i]], method = method)
    jaccard <- getJaccard(sub.kmeans[[i]], method = method)
    fmi <- getFMindex(sub.kmeans[[i]], method = method)
    nmi <- getNMI(sub.kmeans[[i]], method = method)

    cur.dat <- data.frame(
      val = c(ari, jaccard, fmi, nmi),
      # eval = c(rep("ARI", length(ari)), rep("Jaccard", length(jaccard)), rep("FMindex", length(fmi)), rep("NMI", length(nmi))),
      eval <- rep(c("ARI", "Jaccard", "FMindex", "NMI"), each = length(ari)),
      # dist = rep(c("euclidean", "pearson", "spearman"), 4),
      dist = rep(names(ari), 4),
      rep = rep(i, length(c(ari, jaccard, fmi, nmi)))
    )
    # cur.dat <- data.frame(
    #   val = c(ari, jaccard, fmi, nmi),
    #   eval = c(rep("ARI", length(ari)), rep("Jaccard", length(jaccard)), rep("FMindex", length(fmi)), rep("NMI", length(nmi))),
    #   dist = rep(c("euclidean", "manhattan", "correlation", "spearman", "maximum"), 4),
    #   rep = rep(i, length(c(ari, jaccard, fmi, nmi)))
    # )
    result.dat <- rbind(result.dat, cur.dat)
  }
  return(result.dat)
}
