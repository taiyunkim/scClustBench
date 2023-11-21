

######################################################################
## Helper function
######################################################################
"toNumeric" <- function(mat) {
  mat <- apply(mat, 2, as.numeric)
  return (mat)
}


# returns filtered matrix
"filterMatrix" <- function(mat, gene_p = 0.8, cell_p = 0.99) {

  if (class(mat) == "data.frame") {
    mat <- as.matrix(mat)
  }
  if (typeof(mat) == "character") {
    mat <- toNumeric(mat)
  }

  # filter genes
  del <- rowSums(mat == 0) / ncol(mat) > gene_p
  if (any(del == TRUE)) {
    mat <- mat[-which(del == TRUE),]
  }

  # filter cells
  del <- colSums(mat == 0) / nrow(mat) > cell_p
  if (any(del == TRUE)) {
    mat <- mat[,-which(del == TRUE)]
  }
  return (mat)
}

"processMatrix" <- function(mat, count = F) {
  if (count) {
    mat <- cpm(mat)
  }
  mat <- log2(mat + 1)
  return (mat)
}


# return subsets of the matrix
"subsetMatrix" <- function(mat, p = 0.8, rep = 5, seed = 1) {
  # subset the matrix
  tmp <- colnames(mat)
  mat <- mat[,order(tmp)]
  colnames(mat) <- tmp[order(tmp)]
  subsetIndex <- c()
  sub.mat <- list()

  set.seed(seed)
  for (i in 1:rep) {
    subsetIndex[i] <- caret::createDataPartition(colnames(mat), p = p)
    sub.mat[[i]] <- mat[,subsetIndex[[i]]]
  }
  return (sub.mat)
}


"getARI" <- function(obj.list, method = "kmeans") {
  if (method == "kmeans") {
    truth <- obj.list[["truth"]]
    obj.list["truth"] <- NULL
    results <- sapply(obj.list, FUN = function(i) {
      adjustedRandIndex(truth, i$cluster)
    }, USE.NAMES = T)
    # results <- c(
    #   adjustedRandIndex(obj.list$truth, obj.list$euclidean$cluster),
    #   adjustedRandIndex(obj.list$truth, obj.list$manhatta$cluster),
    #   adjustedRandIndex(obj.list$truth, obj.list$correlation$cluster),
    #   adjustedRandIndex(obj.list$truth, obj.list$spearman$cluster),
    #   adjustedRandIndex(obj.list$truth, obj.list$maximum$cluster)
    # )
  } else if (method == "simlr") {
    truth <- obj.list[["truth"]]
    obj.list["truth"] <- NULL
    results <- sapply(obj.list, FUN = function(i) {
      adjustedRandIndex(truth, i$y$cluster)
    }, USE.NAMES = T)

  }

  return (results)
}

"getJaccard" <- function(obj.list, method = "kmeans") {
  truth <- obj.list[["truth"]]
  obj.list["truth"] <- NULL
  if (method == "kmeans") {
    results <- sapply(obj.list, FUN = function(i) {
      cluster_similarity(as.numeric(factor(truth)), as.numeric(factor(i$cluster)), similarity = "jaccard")
    }, USE.NAMES = T)
  } else if (method == "simlr") {
    results <- sapply(obj.list, FUN = function(i) {
      cluster_similarity(as.numeric(factor(truth)), as.numeric(factor(i$y$cluster)), similarity = "jaccard")
    }, USE.NAMES = T)
  }

  return (results)
}

"getFMindex" <- function(obj.list, method = "kmeans") {
  truth <- obj.list[["truth"]]
  obj.list["truth"] <- NULL
  if (method == "kmeans") {
    results <- sapply(obj.list, FUN = function(i) {
      FM_index(truth, i$cluster)
    }, USE.NAMES = T)

  } else if (method == "simlr") {
    results <- sapply(obj.list, FUN = function(i) {
      FM_index(truth, i$y$cluster)
    }, USE.NAMES = T)
  }

  return (results)
}

"getNMI" <- function(obj.list, method = "kmeans") {
  truth <- obj.list[["truth"]]
  obj.list["truth"] <- NULL
  if (method == "kmeans") {
    results <- sapply(obj.list, FUN = function(i) {
      igraph::compare(as.numeric(factor(i$cluster)), as.numeric(factor(truth)), method = "nmi")
    }, USE.NAMES = T)

  } else if (method == "simlr") {
    results <- sapply(obj.list, FUN = function(i) {
      igraph::compare(as.numeric(factor(i$y$cluster)), as.numeric(factor(truth)), method = "nmi")
    }, USE.NAMES = T)
  }

  return (results)
}


# convert a data.frame to matrix where n is the number of repeats
"summarise_dat_to_mat" <- function(dat, n = 5) {
  median_mat <- matrix(nrow = length(unique(dat$dist)), ncol = 1)
  sd_mat <- matrix(nrow = length(unique(dat$dist)), ncol = 1)
  mean_mat <- matrix(nrow = length(unique(dat$dist)), ncol = 1)

  for (j in 1:length(unique(dat$dist))) {
    mat <- dat[dat$dist == unique(dat$dist)[j],]
    median_mat[j,1] <- apply(mat[, 1, drop = F], 2, median)
    sd_mat[j,1] <- apply(mat[, 1, drop = F], 2, sd)
    mean_mat[j,1] <- apply(mat[, 1, drop = F], 2, mean)
  }

  rownames(median_mat) <- unique(dat$dist)
  rownames(sd_mat) <- unique(dat$dist)
  rownames(mean_mat) <- unique(dat$dist)

  se_mat <- sd_mat/sqrt(n)

  return (
    list(
      median = median_mat,
      sd = sd_mat,
      mean = mean_mat,
      se = se_mat
    )
  )
}




######################################################################
## Kmeans
######################################################################


# returns Kmeans results
"getKmeans" <- function(mat, nCs, similarity = NULL, seed = 1, cores = 1, os = "other", ...) {
  # nCluster <- length(levels(factor(colnames(mat))))
  nCluster <- nCs

  if (is.null(similarity)) {
    distMethods <- c("euclidean", "manhattan", "correlation", "spearman", "maximum")
  } else {
    distMethods = similarity
  }

  # kmeans
  if (os == "other") {
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    allResults <- foreach::foreach(
      i = seq_along(distMethods),
      .final = function(x) {
        setNames(x, distMethods)
      }
    ) %dopar%
    {
      set.seed(seed)
      amap::Kmeans(t(mat), nCluster, method = distMethods[i], ...)
    }
    stopCluster(cl)
  } else {
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    allResults <- foreach::foreach(
      i = seq_along(distMethods),
      .final = function(x) {
        setNames(x, distMethods)
      }
    ) %dopar%
    {
      set.seed(seed)
      amap::Kmeans(t(mat), nCluster, method = distMethods[i], ...)
    }
    stopCluster(cl)
  }
  allResults[["truth"]] <- colnames(mat)

  return(allResults)
}

# assuming the colnames is the labels
# subsets the matrix and return their kmeans results
"kmeanSubMatrix" <- function(mat, nCs, p = 0.8, similarity = NULL, rep = 5, seed = 1, cores = 1, os = "other", ...) {
  sub.mat <- subsetMatrix(mat, p, rep, seed)

  allResults <- lapply(sub.mat, FUN = function(i) {
    getKmeans(i, nCs, similarity = similarity, seed = seed, cores = cores, os = os, ...)
  })
  return (allResults)
}


"kmeanSubMatrix_Dist" <- function(mat, p = 0.8, rep = 5, seed = 1, cores = 1,
                                distMethod = c("euclidean", "manhattan", "correlation", "spearman", "maximum")) {
  sub.mat <- subsetMatrix(mat, p, rep, seed)
  results <- mclapply(1:length(sub.mat), mc.cores = cores, function(i) {
    mat <- sub.mat[[i]]
    nCluster <- length(levels(factor(colnames(mat))))
    set.seed(seed)
    result <- Kmeans(t(mat), nCluster, iter.max = 100, nstart = 25, method = distMethod[1])
    list(
      truth = colnames(mat),
      result = result
    )
  })
  return(results)
}



# evaluate cluster results
# assume that the input 'kmeans.obj' is a return object from kmeanSubMatrix
"evaluateKmeans" <- function(sub.kmeans) {
  # for each subset matrix kmeans
  result.dat <- data.frame()
  for (i in 1:length(sub.kmeans)) {
    ari <- getARI(sub.kmeans[[i]])
    jaccard <- getJaccard(sub.kmeans[[i]])
    fmi <- getFMindex(sub.kmeans[[i]])
    nmi <- getNMI(sub.kmeans[[i]])
    cur.dat <- data.frame(
      val = c(ari, jaccard, fmi, nmi),
      eval = c(rep("ARI", length(ari)), rep("Jaccard", length(jaccard)), rep("FMindex", length(fmi)), rep("NMI", length(nmi))),
      dist = rep(c("euclidean", "manhattan", "correlation", "spearman", "maximum"), 4),
      rep = rep(i, length(c(ari, jaccard, fmi, nmi)))
    )
    result.dat <- rbind(result.dat, cur.dat)
  }
  return(result.dat)
}



"run_kmeans_pipeline" <- function(mat, count = F, logMatrix = T, rep = 5, cores = 1) {
  mat <- filterMatrix(mat)
  if (logMatrix) {
    mat <- processMatrix(mat, count = count)
  }
  result <- kmeanSubMatrix(mat)
  return (result)
}





######################################################################
## SIMLR
######################################################################

"getSIMLR" <- function(mat, nCs, similarity = NULL, seed = 1, cores = 1, os = "other", ...) {
  # nCluster <- length(levels(factor(colnames(mat))))
  nCluster <- nCs
  if (is.null(similarity)) {
    distMethods <- c("euclidean", "pearson", "spearman")
  } else {
    distMethods = similarity
  }

  # SIMLR
  if (os == "other") {
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    allResults <- foreach::foreach(
      i = seq_along(distMethods),
      .final = function(x) {
        setNames(x, distMethods)
      }
    ) %dopar%
    {
      set.seed(seed)
      SIMLR(mat, nCluster, distMethod = distMethods[i], ...)
    }
    stopCluster(cl)
  } else {
    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoParallel(cl)
    allResults <- foreach::foreach(
      i = seq_along(distMethods),
      .final = function(x) {
        setNames(x, distMethods)
      }
    ) %dopar%
    {
      set.seed(seed)
      SIMLR(mat, nCluster, distMethod = distMethods[i], ...)
    }
    stopCluster(cl)
  }
  allResults[["truth"]] <- colnames(mat)

  return(allResults)
}



"simlrSubMatrix" <- function(mat, nCs, p = 0.8, similarity = NULL, rep = 5, seed = 1, cores = 1, os = "other", ...) {
  sub.mat <- subsetMatrix(mat, p, rep, seed)

  allResults <- lapply(sub.mat, FUN = function(i) {
    getSIMLR(i, nCs, similarity = similarity, seed = seed, cores = cores, os = os, ...)
  })
  # sub.mat <- subsetMatrix(mat, p = p, rep = rep, seed)
  #
  # if (os == "other") {
  #   cl <- parallel::makeCluster(cores)
  #   doParallel::registerDoParallel(cl)
  #   allResults <- foreach(
  #     i = seq_along(sub.mat)
  #   ) %dopar%
  #   {
  #     getSIMLR(sub.mat[[i]], similarity = similarity, seed = seed, ...)
  #   }
  #   stopCluster(cl)
  # } else {
  #   cl <- parallel::makeCluster(cores)
  #   doSNOW::registerDoSNOW(cl)
  #   allResults <- foreach(
  #     i = seq_along(sub.mat)
  #   ) %dopar%
  #   {
  #     getSIMLR(sub.mat[[i]], similarity = similarity, seed = seed, ...)
  #   }
  #   stopCluster(cl)
  # }

  # results <- getSIMLR(sub.mat, seed = seed)
  return (allResults)
}


