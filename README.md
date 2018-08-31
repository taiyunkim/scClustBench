# scClust

This package is to perform clustering on single-cell RNA-sequencing data with different similarity metrics on k-means and SIMLR clustering methods.

## Requirements

Before the installation, some of R packages must be installed.
These include:

* parallel
* amap
* caret
* mclust
* clusteval
* dendextend
* igraph
* foreach
* doParallel
* doSNOW
* ggplot2
* Hmisc 


## Installation


```r
devtools::install_github("taiyunkim/scClustBench", build_vignettes = TRUE)
library(scClust)
```

Building the vignette may take some time. If you wish not to create the vignette during installation, try:

```r
devtools::install_github("taiyunkim/scClustBench")
library(scClust)
```

**NOTE:** *For mac users, the official cran mirror of R tools for OS X and R tools for OS X on r.research.att.com that lists the gfortran binary are out of date. You will need to update `gfortran` prior to this package installation.* 

Current version of this package is implemented to run SIMLR (Wang et al, 2017) or k-means clustering methods with various similarity metrics.

Available metrics include:

&nbsp;&nbsp;&nbsp;&nbsp;**SIMLR** - `"pearson"` correlation, `"spearman"` correlation and `"euclidean"` distances.

&nbsp;&nbsp;&nbsp;&nbsp;**K-means** - `"pearson"` correlation, `"spearman"` correlation, `"euclidean"` distance, `"manhattan"` distance and `"maximum"` distance.



## Using scClust

For further demonstration, see:

```r
browseVignette("scClust")
```

### Section 1. Clustering with different similarity metrics with `scClust`

To run scClust, 

```r
# SIMLR
simlr.result <- scClust(mat, nCs, similarity = "pearson", method = "simlr", seed = 1, cores.ratio = 0)

# K-means
kmeans.result <- scClust(mat, nCs, similarity = "pearson", method = "kmeans", seed = 1, nstart = 10, iter.max = 10)
```

This function allows you to perform clustering with a user specified similarity metrics. The return values of `scClust` are identical to clustering methods for `Kmeans` and `SIMLR` functions.


### Section 2. Benchmarking different similarity metrics with `scClustBench`

This section is to compare a set of similarity metrics on clustering methods to benchmark their perfomance accuracy.

To run scClustBench,

```r
# SIMLR
simlr.result <- scClustBench(mat, method = "simlr", rep = 2, seed = 1, cores = 1, cores.ratio = 0)

# K-means
kmeans.result <- scClustBench(mat, method = "kmeans", rep = 2, seed = 1, cores = 1, nstart = 10, iter.max = 10)
```


You can evaluate this result with the function `evalScClustBench` and plot with `plotSimlrEval` or `plotKmeansEval`.




