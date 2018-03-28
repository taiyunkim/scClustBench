# Impacts of similarity metrics on scRNA-seq

## Overview
A common goal in scRNA-seq data analysis is to discover and characterise cell types typically through clustering. The quality of the clustering output therefore plays a critical role in biological discovery. While numerous clustering algorithms were proposed for scRNA-seq data analysis, fundamentally all clustering algorithms rely on a similarity metric for categorising individual cells. 

Here, we have utilised k-means clustering from <a href="https://cran.r-project.org/web/packages/amap/index.html">amap</a> R package and also implemented the state-of-the-art kernel-based clustering algorithm <a href="https://bioconductor.org/packages/release/bioc/html/SIMLR.html">SIMLR</a> to investigate the effect of different similarity metrics in clustering scRNA-seq results.


## R Packages required

<ul>
	<li>mySIMLR (installation required using 'mySIMLR_1.4.0.tar.gz' file)</li>
	<li>caret</li>
	<li>parallel</li>
	<li>edgeR</li>
	<li>mclust</li>
	<li>clusteval</li>
	<li>dendextend</li>
	<li>NMF</li>
	<li>igraph</li>
</ul>

## mySIMLR package installation

### In Terminal,

1. Open terminal
2. Type the following line:
	
	```bash
	$ R CMD INSTALL /PATH/TO/FILE/mySIMLR_1.4.0.tar.gz
	```
	
### In R,
1. Start R
2. Type the following line in console:

	```r
	> install.packages("/PATH/TO/FILE/mySIMLR_1.4.0.tar.gz", repos = NULL, type = "source")
	```


## Examples

### Load Data and functions

```r
source("functions.R")
load("singleCellData.RData")
```
Load your single cell data say `singleCellData.RData` for this example to the environment. This is m*n matrix where n is the number of cells and m is the number of genes.

### Data Preprocessing
`singleCellData` is an expression matrix that contains a lot of zeros in the data and hence require processing to remove genes that are very lowly expressed. The quantification of the expression is in CPM (counts per million) and requires log before clustering.

```r
mat <- singleCellData
# Filtering step
mat.filtered <- filterMatrix(mat)

# log matrix
mat.filtered.log <- logMatrix(mat.filtered, count = F)
```


### Clustering

Perform clustering on processed data using Kmeans or SIMLR.

* K-means
	
	```r
	kmeans.result <- kmeanSubMatrix(mat.filtered.log, cores = 1)
	```

* SIMLR

	```r
	simlr.result <- simlrSubMatrix(mat.filtered.log, cores.ratio = 1)
	```

### Evaluate and plot

To evaluate and visualise your results from clustering, use the following functions for each methods.

* K-means

	```r
	# evaluate kmeans result data.frame
	kmeans.eval <- evaluateKmeans(kmeans.result)
	
	# plot evaluation
	p <- plotKmeansEval(kmeans.eval)

	```

* SIMLR

	```r
	# evaluate SIMLR result data.frame
	simlr.eval <- evaluateSIMLR(simlr.result)

	# plot evaluation
	p <- plotSimlrEval(simlr.eval)
	```



## Contacts
* <a href=mailto:taiyun.kim@sydney.edu.au>Taiyun Kim</a>
* <a href=mailto:pengyi.yang@sydney.edu.au>Dr. Pengyi Yang</a>
* <a href=mailto:jean.yang@sydney.edu.au>Prof. Jean Yang</a>


