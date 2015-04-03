\name{TAPPA}
\alias{TAPPA}

\title{
Function to use TAPPA method on microarray or RNA-Seq data
}
\description{
The functions analyses the differential expression of pathways via TAPPA method. Expression is compared between two groups of samples by Mann-Whitney test. P-values are later adjusted for multiple hypothesis testing by Benjamini-Hochberg's FDR method.
}
\usage{
TAPPA(x, group, pathways, type, preparePaths=TRUE, norm.method=NULL, test.method=NULL, test=t.test, normalize=TRUE, verbose=FALSE, both.directions=TRUE, 
   maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="entrez", convertBy=NULL)
}

\arguments{
  \item{x}{
An \code{ExpressionSet} object or a gene expression data matrix or count matrix, rows refer to genes, columns to samples
}
  \item{group}{
Name or number of the phenoData column or a character vector or factor that contains required class assigments
}
  \item{pathways}{
A list of pathways in a form from \code{graphite} package or created by \code{preparePathways()}
}
  \item{type}{
Type of the input data, \code{"MA"} for microarray and \code{"RNASeq"} for RNA-Seq
}
\item{preparePaths}{
Logical, by default the pathways are transformed with \code{preparePathways()}. Use \code{FALSE}, if you have done this transformation separately}
\item{norm.method}{Character, the method to normalize RNAseq data. If \code{NULL} then TMM-normalization is performed. Possible values are: \code{"TMM", "DESeq2", "rLog", "none"}}
\item{test.method}{Character, the method for differentiall expression analysis of RNAseq data. If \code{NULL} then \code{"voomlimma"} is used. Possible values are: \code{"DESeq2", "voomlimma", "vstlimma", "edgeR"}. This analysis is needed only for the visualization.   }
\item{test}{Function implementing a statistical test comparing PCI scores between groups. It is employed as \code{test(PCI~group)$p.value}, where \code{PCI} is a numeric vector of the same length as \code{group}}
  \item{normalize}{
Logical, should data be normalized? 
}
  \item{verbose}{
Logical, if \code{TRUE} names of the pathways are printed as they are analysed
}
  \item{both.directions, maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy}{Arguments for the \code{preparePathways()}}
}

\value{ 
A list,
\item{res}{A data frame, rows refer to pathways. Columns contain: number of valid PCI-scores, median, min and max of the PCI scores for each group of samples, p-value of the \code{test} (\code{p.val}) and adjusted p-value (\code{p.adj}). If less than two nodes are present in the data, the function puts \code{NA}'s in all columns. }
\item{topo.sig}{\code{NULL}, it is preserved for the compatibility with other methods implemented in this package}
\item{degtest}{A numeric vector of gene-level differential expression statistics}
}
\references{
Gao, S. and Wang, X. (2007) TAPPA: topological analysis of pathway phenotype association. Bioinformatics, 23, pages 3100-3102
}
\author{Ivana Ihnatova}
\seealso{
 \code{\link{preparePathways}}}
\examples{

if (require(DEGraph)) {
  data("Loi2008_DEGraphVignette")
  TAPPA(exprLoi2008, classLoi2008, biocarta[1:10], type="MA")
}

\dontrun{
if (require(gageData)) {

 data(hnrnp.cnts)
 group<-c(rep("sample",4), rep("control",4))
 hnrnp.cnts<-hnrnp.cnts[rowSums(hnrnp.cnts)>0,]
 TAPPA( hnrnp.cnts, group, biocarta[1:10], type="RNASeq", norm.method="TMM")
}
}
}


\keyword{ htest }

