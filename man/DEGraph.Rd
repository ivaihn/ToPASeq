\name{DEGraph}
\alias{DEGraph}

\title{
Function to use DEGraph method on microarray or RNA-Seq data
}
\description{
DEGraph implements recent hypothesis testing methods which directly assess whether a particular gene network is differentially expressed between two conditions. In employs Graph Laplacian, Fourier transformation and multivariate T2-statistic
}
\usage{
DEGraph(x, group, pathways, type, preparePaths=TRUE, norm.method=NULL, test.method=NULL, overall="biggest", useInteractionSigns=TRUE, EdgeAttrs=NULL,
   both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="entrez", convertBy=NULL)
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
Type of the data, \code{"MA"} for microarray and \code{"RNASeq"} for RNA-Seq
}
\item{preparePaths}{
Logical, by default the pathways are transformed with \code{preparePathways()}. Use \code{FALSE}, if you have done this transformation separately}
\item{norm.method}{Character, the method to normalize RNAseq data. If \code{NULL} then TMM-normalization is performed. Possible values are: \code{"TMM", "DESeq2", "rLog", "none"}}
\item{test.method}{Character, the method for differentiall expression analysis of RNAseq data. If \code{NULL} then \code{"voomlimma"} is used. Possible values are: \code{"DESeq2", "voomlimma", "vstlimma", "edgeR"}. This analysis is needed only for the visualization.   }
  \item{overall}{Character, how should the overall p-value for a pathway be calculated. The possible values are: "mean", "min", "biggest". "biggest" returns the p-value of the biggest connected component.  }
 \item{useInteractionSigns}{Logical, should types of interaction be included in the analysis?}
\item{EdgeAttrs}{A list containing two data.frames. See \code{makeDefaultEdgeData()} for the details. The interactions are assigned signs according to the \code{beta} column of the second data.frame. The procedure is similar to the SPIA method}
\item{both.directions, maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy}{Arguments for the \code{preparePathways()}}
}

\value{
A list:
\item{res}{Results from analysis of individual pathways. The first column refers to the overall p-value for a pathway. Then groups of four columns follows. One group refers to one connected component and contains a pair of p-values (without and with Fourier transformation), graph and number of Fourier componets used in the test. The number of groups is equal to the highest number of components in analysed pathways. Components are sorted in the decreasing order of their nodes number.}
 \item{topo.sig }{\code{NULL}, present for the compatibility with outputs from other methods}
 \item{degtest}{A data.frame of gene-level statistics of all genes in the dataset}
}
\references{
L. Jacob, P. Neuvial, and S. Dudoit. Gains in power from structured two-sample tests of means on graphs. Technical Report arXiv:q-bio/1009.5173v1, arXiv, 2010. 
}
\author{
Ivana Ihnatova
}
\seealso{
 \code{\link{preparePathways}}}
\examples{

if (require(DEGraph)) {
  data("Loi2008_DEGraphVignette")

  DEGraph(exprLoi2008, classLoi2008, biocarta[1:10], type="MA")
}
\dontrun{
if (require(gageData)) {

 data(hnrnp.cnts)
 hnrnp.cnts<-hnrnp.cnts[rowSums(hnrnp.cnts)>0,]
 group<-c(rep("sample",4), rep("control",4))
 DEGraph(hnrnp.cnts, group, biocarta[1:10], type="RNASeq", norm.method="TMM")
}
}
}
\keyword{ htest }

