# QDNAseq.hg38
>QDNAseq bin annotation for hg38

[![HitCount](http://hits.dwyl.io/asntech/QDNAseq.hg38.svg)](http://hits.dwyl.io/asntech/QDNAseq.hg38)

This package is a fork of Bioconductor R package
[QDNAseq.hg19](https://doi.org/doi:10.18129/B9.bioc.QDNAseq.hg19)
This package provides QDNAseq bin annotations of size `5, 10, 15, 30, 50, 100, 500 and 1000` kbp for the human genome build hg38

## Installation

Install `hg38` package from GitHub:

``` r
#Install the QDNAseq.hg38 package using remotes
remotes::install_github("asntech/QDNAseq.hg38@main")
#or devtools
devtools::install_github("asntech/QDNAseq.hg38@main")
```

## Use QDNAseq.hg38

``` r
library(QDNAseq)
library(QDNAseq.hg38)
bins <- getBinAnnotations(binSize=50, genome="hg38")
```

Find more details about QDNAseq here: https://doi.org/doi:10.18129/B9.bioc.QDNAseq