# QDNAseq.hg38
>QDNAseq bin annotation for hg38

This package provides QDNAseq bin annotations of size `1, 5, 10, 15, 30, 50, 100, 500 and 1000` kbp for the human genome build hg38.The bin annotations are created using the steps mentioned in QDNAseq vignette and also [here](https://github.com/ccagc/QDNAseq/issues/59).


## Installation

Install the package from GitHub:

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

`QDNAseq.hg38` is adopted from [QDNAseq.hg19](https://doi.org/doi:10.18129/B9.bioc.QDNAseq.hg19). Find more details about QDNAseq here: https://doi.org/doi:10.18129/B9.bioc.QDNAseq
