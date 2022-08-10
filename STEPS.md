## Steps used to generate hg38 bins

The following steps were used to create the bin files.

### 1. Create mappability file
To calculate the average mappabilities, we need a mappability file in the `bigWig` format and the `bigWigAverageOverBed` binary.

We used [GenMap](https://github.com/cpockrandt/genmap) to generate the 50mer mappability file with 2-mismatches.

``` bash
# Download the pre-build index for GRCh38/hg38
wget http://ftp.imp.fu-berlin.de/pub/cpockrandt/genmap/indices/grch38-no-alt.tar.gz

# Compute 50mer mappability file in wig format
genmap map -K 50 -E 2 -I /path/to/index/grch38-no-alt -O /path/to/output/folder -w

# Convert the wig file to bigwig
wigToBigWig grch38-no-alt.wig <hg38.chrom.sizes> mappability.genmap.50mer.bigwig

```

### 2. Get ENCODE excluded regions

``` bash
# Download ENCODE excluded regions aka blacklist regions
wget https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz?raw=true -o hg38-blacklist.v2.bed.gz
gzip -d hg38-blacklist.v2.bed.gz
```

### 3. Control data set to calculate median residuals

To calculate median residuals of the LOESS fit a control set of 38 samples from the 1000 Genomes Project were downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/.

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

Next, reads were trimmed to 50 bp using `trim_galore v0.4.4`, and the multiple files for each sample.

``` bash
for read in *.fq.gz; do
  trim_galore --cores 32 --hardtrim5 50 $read
done;
```

Next, these reads were aligned with `BWA v0.6.2` allowing two mismatches and end-trimming of bases
with qualities below 40. bam files were sorted and indexed using `samtools v1.11`

``` bash
$REF_GENOME=/path/to/bwa/idex/
for SAMPLE in *.fq.gz; do
  bwa aln -t 16 -n 2 -q 40 $REF_GENOME ${SAMPLE}.fq.gz > ${SAMPLE}.sai
  bwa samse $REF_GENOME ${SAMPLE}.sai ${SAMPLE}.fq.gz | samtools view --threads 16 -hb | samtools sort - > ${SAMPLE}.bam
  samtools index ${SAMPLE}.bam
  rm ${SAMPLE}.sai
done;
```

### Create the bin annotations

Next we used `R v3.6.0` with `BSgenome.Hsapiens.UCSC.hg38` and `QDNAseq` to generate bins of size `1, 5, 10, 15, 30, 50, 100, 500 and 1000` kbp.

``` r

library(Biobase)
library(BSgenome.Hsapiens.UCSC.hg38)
library(QDNAseq)
library(future)

#set virtual mem
options(future.globals.maxSize= 8912896000)

# change the current plan to access parallelization
future::plan("multiprocess", workers = 4)

for (binsize in c(1000, 500, 30, 15, 50, 10, 5, 1)) {

  bins <- createBins(bsgenome=BSgenome.Hsapiens.UCSC.hg38, binSize=binsize)
  bins$mappability <- calculateMappability(bins,
    bigWigFile="/path/to/hg38/mappability.genmap.50mer.bigwig",
    bigWigAverageOverBed="/path/to/bigWigAverageOverBed")

  bins$blacklist <- calculateBlacklist(bins, bedFiles=c(
    "/path/to/hg38-blacklist.v2.bed"))

  bins$residual <- NA
  bins$use <- bins$chromosome %in% as.character(1:22) & bins$bases > 0
  
  #
  tg <- binReadCounts(bins,
    path="/path/to/1000Genomes/hg38/bams", cache=TRUE)

  bins$residual <- iterateResiduals(tg)
  
  bins <- AnnotatedDataFrame(bins,
    varMetadata=data.frame(labelDescription=c(
    "Chromosome name",
    "Base pair start position",
    "Base pair end position",
    "Percentage of non-N nucleotides (of full bin size)",
    "Percentage of C and G nucleotides (of non-N nucleotides)",
    "Average mappability of 50mers with a maximum of 2 mismatches",
    "Percent overlap with ENCODE blacklisted regions",
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
    "https://github.com/asntech/QDNAseq.hg38/raw/master/data/hg38.",
    binsize, "kbp.SR50.rda"),
    md5=digest::digest(bins@data),
    sessionInfo=sessionInfo())

  attr(bins, "QDNAseq") <- QDNAseqInfo
  save(bins, file=paste0("hg38.", binsize, "kbp.SR50.rda"), compress='xz')
}

```