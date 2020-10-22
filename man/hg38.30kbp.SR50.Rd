\name{hg38.30kbp.SR50}
\alias{hg38.30kbp.SR50}
\docType{data}

\title{hg38 30kbp bin annotations}

\description{ Bin annotations are caclulated for non overlapping 30kbp bins 
generated as described in Scheinin et al. (see references). The annotated data 
frame contains:
\itemize{
\item chromosome: Chromosome name,
\item start: Base pair start position,
\item end: Base pair end position,
\item bases: Percentage of non-N nucleotides (of full bin size),
\item gc: Percentage of C and G nucleotides (of non-N nucleotides),
\item mappability: Average mappability of 50mers with a maximum of 2 
mismatches as described in by Derrien et al. (see references),
\item blacklist: Percent overlap with ENCODE blacklisted regions (see 
references), 
\item residual: Median loess residual calculated from 1000 Genomes (see 
references), 
\item use: Whether the bin should be used in subsequent analysis steps,
}
}

\examples{
data("hg38.30kbp.SR50")
assign("bins", get("hg38.30kbp.SR50"))
\dontrun{readCounts <- binReadCounts(bins=bins, path="./bam")}

# or

bins <- getBinAnnotations(binSize=30, genome="hg38")
\dontrun{readCounts <- binReadCounts(bins=bins, path="./bam")}
}

\value{ Returns an \code{\link[Biobase]{AnnotatedDataFrame}} object. }

\author{Daoud Sie}
\keyword{datasets}

\references{
DNA copy number analysis of fresh and formalin-fixed specimens by shallow 
whole-genome sequencing with identification and exclusion of problematic 
regions in the genome assembly. Scheinin I, Sie D, Bengtsson H, van de Wiel M, 
Olshen A, van Thuijl H, van Essen H, Eijk P, Rustenburg F, Meijer G, Reijneveld
J, Wesseling P, Pinkel D, Albertson D, Ylstra B 2014 \emph{Genome Research}
vol: \bold{24 (12)} pp: 1--11

Fast Computation and Applications of Genome Mappability. Derrien T, Estelle J,
Sola S, Knowles D, Raineri E, Guigo R, Ribeca P January 19, 2012
\emph{PLOS ONE} doi: 10.1371/journal.pone.0030377

An integrated map of genetic variation from 1,092 human genomes.
1000 Genomes Project Consortium, Abecasis GR, Auton A, Brooks LD, DePristo MA,
Durbin RM, Handsaker RE, Kang HM, Marth GT, McVean GA 2012
\emph{Nature} Nov 1; \bold{491(7422)}:56--65.

An integrated encyclopedia of DNA elements in the human genome.
ENCODE Project Consortium 2012
\emph{Nature} Sep 6; \bold{489(7414)}:57--74.
}
