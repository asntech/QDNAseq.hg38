\name{QDNAseq.hg38-package}
\alias{QDNAseq.hg38-package}
\alias{QDNAseq.hg38}
\docType{package}
\title{Package QDNAseq.hg38}
\description{
This package provides QDNAseq binannotations for the mouse genome build hg38
for bin sizes 1, 5, 10, 15, 30, 50, 100, 500 and 1000 kbp (kilobasepair).

The datasets are named as follows:

hg38.1kbp.SR50\cr
hg38.5kbp.SR50\cr
hg38.10kbp.SR50\cr
hg38.15kbp.SR50\cr
hg38.30kbp.SR50\cr
hg38.50kbp.SR50\cr
hg38.100kbp.SR50\cr
hg38.500kbp.SR50\cr
hg38.1000kbp.SR50\cr
}

\examples{
data("hg38.50kbp.SR50")
assign("bins", get("hg38.50kbp.SR50"))
\dontrun{readCounts <- binReadCounts(bins=bins, path="./bam")}

# or

bins <- getBinAnnotations(binSize=50, genome="hg38")
\dontrun{readCounts <- binReadCounts(bins=bins, path="./bam")}
}

\author{Daoud Sie}
\section{License}{
This package is licensed under GPL.
}
\keyword{package}

