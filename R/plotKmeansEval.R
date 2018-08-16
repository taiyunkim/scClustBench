#' plot evaluations returned from \emph{evalScClustBench}
#'
#' @title plotKmeansEval
#'
#' @examples
#'
#'
#' @param kmeans.eval An object returned from \emph{evalScClustBench}
#' @param eval.method Evaluation method to be plotted. Available options are: NMI (default), FM, Jaccard and ARI.
#'
#'
#' @return ggplot object
#'
#' @export plotKmeansEval
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
#' @import ggplot2
#' @useDynLib scClust projsplx
#'
"plotKmeansEval" <- function(kmeans.eval, eval.method = c("NMI", "FM", "ARI", "Jaccard")) {
  dat <- kmeans.eval[kmeans.eval$eval == eval.method[1],]
  dat.summary <- summarise_dat_to_mat(dat)

  med <- reshape2::melt(dat.summary$median)
  se <- reshape2::melt(dat.summary$se)
  med <- cbind(med,se)

  colnames(med) <- c("dist1", "data1", "median", "dist2", "data2", "se")

  med <- cbind(med, rep(NA, nrow(med)))
  colnames(med) <- c(colnames(med)[-ncol(med)], "type")

  tmp <- med

  med$type[(med$dist1 == "correlation") | (med$dist1 == "spearman")] <- "correlation"
  med$type[!((med$dist1 == "correlation") | (med$dist1 == "spearman"))] <- "abs magnitude"


  med_dist <- med[med$type == "abs magnitude",]
  med_shape <- med[med$type == "correlation",]


  levels(med$dist1)[3] = "pearson"

  p <- ggplot(med, aes( x = data1, y = median, fill = dist1)) +
    geom_bar(stat = "identity", position = "dodge")  +
    geom_errorbar(aes(ymin = median-se, ymax = median+se),
                  width = .2,
                  position = position_dodge(.9)) +
    scale_fill_manual("legend",
                      values = c("maximum" = "#feb24c",
                                 "manhattan" = "#fd8d3c", "euclidean" = "#fc4e2a",
                                 "spearman" = "#a6bddb", "pearson" = "#054287")) +
    labs(title = "NMI bar plot (median + se)", x = "Datasets", y = "NMI values") +
    scale_y_continuous(limits = c(0, 1)) +
    coord_flip()
  return (p)

}
