#' Evaluate cluster result from \emph{scClustBench}
#'
#' @title evaluateSIMLR
#'
#' @examples
#' evaluateSIMLR(sub.simlr)
#'
#' @param sub.simlr A object returned from \emph{scClustBench}
#'
#' @return
#'
#'
#' @export evaluateSIMLR
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
"evaluateSIMLR" <- function(sub.simlr) {
  result.dat <- data.frame()
  for (i in 1:length(sub.simlr)) {
    ari <- getARI(sub.simlr[[i]], method = "simlr")
    jaccard <- getJaccard(sub.simlr[[i]], method = "simlr")
    fmi <- getFMindex(sub.simlr[[i]], method = "simlr")
    nmi <- getNMI(sub.simlr[[i]], method = "simlr")
    cur.dat <- data.frame(
      val = c(ari, jaccard, fmi, nmi),
      eval = c(rep("ARI", length(ari)), rep("Jaccard", length(jaccard)), rep("FMindex", length(fmi)), rep("NMI", length(nmi))),
      dist = rep(c("euclidean", "pearson", "spearman"), 4),
      rep = rep(i, length(c(ari, jaccard, fmi, nmi)))
    )
    result.dat <- rbind(result.dat, cur.dat)
  }
  return(result.dat)
}
