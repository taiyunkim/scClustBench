# scClust

**scClust** is a package which performs k-means and SIMLR clustering on single-cell RNA-sequencing data using various similarity metrics.

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
devtools::install_github("SydneyBioX/scClustBench", build_vignettes = TRUE)
library(scClust)
```

Building the vignette may take some time. If you wish not to create the vignette during installation, try:

```r
devtools::install_github("SydneyBioX/scClustBench")
library(scClust)
```

**NOTE:** *For mac users, the official cran mirror of R tools for OS X and R tools for OS X on r.research.att.com that lists the gfortran binary are out of date. You will need to update `gfortran` and add the following line `FLIBS=-L/usr/local/Cellar/gcc/X.Y.Z/lib/gcc/X` (where `X.Y.Z` is your gcc version) to `~/.R/Makevars` prior to this package installation.* 

Current version of this package is implemented to run SIMLR (Wang et al, 2017) or k-means clustering methods with various similarity metrics.

Available metrics include:

&nbsp;&nbsp;&nbsp;&nbsp;**SIMLR** - `"pearson"` correlation, `"spearman"` correlation and `"euclidean"` distances.

&nbsp;&nbsp;&nbsp;&nbsp;**K-means** - `"pearson"` correlation, `"spearman"` correlation, `"euclidean"` distance, `"manhattan"` distance and `"maximum"` distance.



## Using scClust

For further demonstration, see:

```r
browseVignettes("scClust")
```

### Load Data

```r
data(GSE82187.sample)
mat <- GSE82187

mat <- log2(mat+1)

# set number of clusters (classes defined in colnames)
nCs <- length(table(colnames(mat))
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
simlr.result <- scClustBench(mat, nCs, method = "simlr", rep = 2, seed = 1, cores = 1, cores.ratio = 0)

# K-means
kmeans.result <- scClustBench(mat, nCs, method = "kmeans", rep = 2, seed = 1, cores = 1, nstart = 10, iter.max = 10)
```


You can evaluate this result with the function `evalScClustBench` and plot with `plotSimlrEval` or `plotKmeansEval`.




