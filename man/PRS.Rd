\name{PRS}
\alias{PRS}

\title{
Function to use PRS method on microarray or RNA-Seq data
}
\description{
A function runs PRS method on a gene expression data matrix or count matrix and vector dividing samples into two groups and a set of pathways from \code{graphite} package. The PRS method (please see Reference for the details) was adapted to \code{graphite}'s graphs where each node is represented only by one gene.
}
\usage{
PRS(x, group, pathways, type, preparePaths=TRUE, norm.method=NULL, test.method=NULL, p.th=0.05, logFC.th=2, nperm=1000,
 both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL)
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
Type of the data, \code{"MA"} for microarray, \code{"RNASeq"} for RNA-Seq, \code{DEtable} data.frame from differential expression analysis, or \code{DEGlist} a list of: log fold-changes of differentially expressed genes and names of the all genes analyses
}
\item{preparePaths}{
Logical, by default the pathways are transformed with \code{preparePathways()}. Use \code{FALSE}, if you have done this transformation separately}

\item{norm.method}{Character, the method to normalize RNAseq data. If \code{NULL} then TMM-normalization is performed. Possible values are: \code{"TMM", "DESeq2", "rLog", "none"}. Ignored for type: "MA","DEtable", "DElist"}
\item{test.method}{Character, the method for differentiall expression analysis of RNAseq data. If \code{NULL} then \code{"voomlimma"} is used. Possible values are: \code{"DESeq2", "voomlimma", "vstlimma", "edgeR"}. Ignored for type: "MA","DEtable", "DElist"   }
\item{p.th}{Numeric, threshold for p-values of tests for differential expression of genes. Use \code{1} if you don't want any threshold to be applied}
\item{logFC.th}{Numeric, threshold for log fold-change of a gene to identify the gene as differentially expressed. Use negative if you don't want any threshold to be applied}
\item{nperm}{Numeric, number of permutations}


\item{both.directions, maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy}{Arguments for the \code{preparePathways()}}
}
\value{
A list,
\item{res}{A data frame with normalized score, p-value and FDR-adjusted p-value for each pathway}
\item{topo.sig}{A list with log fold-changes and number of downstream differentially expressed nodes for nodes of individual pathways}
\item{degtest}{A named vector of statistics from testing the differential expression of genes}

}
\references{
Maysson Al-Haj Ibrahim, Sabah Jassim, Michael Anthony Cawthorne, and Kenneth Langlands. A Topology-Based Score for Pathway Enrichment,
 Journal of Computational Biology. May 2012, 19(5): 563-573
}
\author{
Ivana Ihnatova
}
\seealso{
 \code{\link{preparePathways}}}
 
\examples{

if (require(DEGraph)) {
  data("Loi2008_DEGraphVignette")
  pathways<-pathways("hsapiens","biocarta")[1:10]
  PRS( exprLoi2008, classLoi2008, pathways, type="MA",  logFC.th=-1, nperm=100)
}
\dontrun{
if (require(gageData)) {

 data(hnrnp.cnts)
 hnrnp.cnts<-hnrnp.cnts[rowSums(hnrnp.cnts)>0,]
 group<-c(rep("sample",4), rep("control",4))
 pathways<-pathways("hsapiens","biocarta")[1:10]
 PRS(hnrnp.cnts, group, pathways, type="RNASeq", logFC.th=-1, nperm=100, test="vstlimma")
}
}
}

\keyword{ htest }

