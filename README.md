[![DOI](https://zenodo.org/badge/306393966.svg)](https://zenodo.org/badge/latestdoi/306393966)

# QDNAseq.hg38: QDNAseq bin annotation for hg38

>QDNAseq bin annotation for the human genome build hg38

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


## Reproducibility

Below are the steps and data used to generate the QDNAseq.hg38 bin annotations.

### Step 1: Baseline samples and processing 

The the fastq files for the following 38 samples was downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/.

These 38 samples are matched the number mentioned in the QDNAseq instructions for hg19 annotations.

``` bash
HG01101
HG01204
HG01495
HG01522
NA06994
NA07051
NA11830
NA11832
NA11918
NA11919
NA11992
NA11994
NA12005
NA18489
NA18498
NA18504
NA18623
NA18624
NA18632
NA18633
NA18636
NA18867
NA18912
NA18960
NA18982
NA18984
NA18986
NA19058
NA19063
NA19064
NA19066
NA19116
NA19138
NA19474
NA19703
NA19707
NA19789
NA19901
```

Next reads for these 38 samples were trimed to 50bp using `TrimGalore v0.4.4`.

```
trim_galore --cores 32 --hardtrim5 50 $FASTQ
```

Next reads were aligned using bwa aln and samse (`BWA v0.6.2`), and sorted and indexed using `Samtools v1.11`

Here is example for one sample `HG01101`.
``` bash
SAMPLE=HG01101
REF_GENOME=/path/to/hg38/genome.fa

bwa aln -t 16 -n 2 -q 40 $REF_GENOME ${SAMPLE}.fq.gz > ${SAMPLE}.sai
bwa samse  ${SAMPLE}.sai ${SAMPLE}.fq.gz | samtools view --threads 16 -hb | samtools sort - > ${SAMPLE}.bam
samtools index ${SAMPLE}.bam
```


### Step 2: Getting mapability bigWig

A 50mer genome mappability file was created using the GenMap (https://github.com/cpockrandt/genmap).
We used the precomputed index file for hg38/GRCh38 was downloaded from `http://ftp.imp.fu-berlin.de/pub/cpockrandt/genmap/indices/grch38-no-alt.tar.gz`

### Step 3: Getting the bin annotations

The following Rscript (`R v3.6`) was used to get the final bins of different sizes (kb). 

``` r
library(Biobase)
library(BSgenome.Hsapiens.UCSC.hg38)
library(QDNAseq)
library(future)

options(future.globals.maxSize= 8912896000)

# change the current plan to access parallelization
future::plan("multiprocess", workers = 4)

for (binsize in c(1000, 500, 30, 15, 50, 10, 5, 1)) {
#for (binsize in c(200,300,400)) {

  bins <- createBins(bsgenome=BSgenome.Hsapiens.UCSC.hg38, binSize=binsize)
  bins$mappability <- calculateMappability(bins,
    bigWigFile="/path/for/bigwig/grch38-no-alt.genmap.50mer.bigwig",
    bigWigAverageOverBed="/path/for/bin/bigWigAverageOverBed")

  bins$blacklist <- calculateBlacklist(bins, bedFiles=c(
    "/path/for/wgEncodeDacMapabilityConsensusExcludable.hg38.bed",
    "/path/for/hg38-blacklist.v2.bed"
    ))

  bins$residual <- NA
  bins$use <- bins$chromosome %in% as.character(1:22) & bins$bases > 0
  
  tg <- binReadCounts(bins,
    path="/path/for/hg38/bams/", cache=TRUE)

  bins$residual <- iterateResiduals(tg)
  
  bins <- AnnotatedDataFrame(bins,
    varMetadata=data.frame(labelDescription=c(
    "Chromosome name",
    "Base pair start position",
    "Base pair end position",
    "Percentage of non-N nucleotides (of full bin size)",
    "Percentage of C and G nucleotides (of non-N nucleotides)",
    "Average mappability of 50mers with a maximum of 2 mismatches",
    "Percent overlap with ENCODE blocklisted regions",
    "Median loess residual from 1000 Genomes (50mers)",
    "Whether the bin should be used in subsequent analysis steps"),
    row.names=colnames(bins)))

  QDNAseqInfo <- list(
    author="Aziz Khan",
    date=Sys.time(),
    organism='Hsapiens',
    build='hg38',
    version=packageVersion("QDNAseq"),
    url=paste0(
    "https://github.com/asntech/QDNAseq.hg38/raw/main/data/hg38.",
    binsize, "kbp.SR50.rda"),
    md5=digest::digest(bins@data),
    sessionInfo=sessionInfo())

  attr(bins, "QDNAseq") <- QDNAseqInfo
  save(bins, file=paste0("hg38.", binsize, "kbp.SR50.rda"))
}

```
