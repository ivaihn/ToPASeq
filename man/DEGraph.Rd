\name{DEGraph}
\alias{DEGraph}

\title{
Function to use DEGraph method on microarray or RNA-Seq data
}
\description{
DEGraph implements recent hypothesis testing methods which directly assess whether a particular gene network is differentially expressed between two conditions. In employs Graph Laplacian, Fourier transformation and multivariate T2-statistic
}
\usage{
DEGraph(x, group, pathways, type = "MA", convert = TRUE, IDs = "entrez", 
        gene.stat="logFC",both.directions=TRUE, overall="biggest",  useInteractionSigns=TRUE, norm.method=NULL)
}
\arguments{
  \item{x}{
An \code{ExpressionSet} object or a gene expression data matrix or count matrix, rows refer to genes, columns to samples
}
  \item{group}{
Name or number of the phenoData column or a character vector or factor that contains required class assigments
}
  \item{pathways}{
A list of pathways in a form from \code{graphite} package
}
  \item{type}{
Type of the data, \code{"MA"} for microarray and \code{"RNASeq"} for RNA-Seq
}
 \item{convert}{
Logical, should nodes' labels be converted?
}
  \item{IDs}{
Character, which gene IDs should be used. One of \code{"entrez"} or \code{"symbol"}
}
  \item{gene.stat}{Character, which gene-level statistic should be exported? \code{"logFC"} stands for log fold-change and \code{"stats"} for the modified t-statistic from \code{limma} package}
 \item{both.directions}{Logical, should undirected edges be substitued by two directed edges (one in each direction) or only by one directed edge oriented as in \code{edges()}}
  \item{overall}{Character, how should the overall p-value for a pathway be calculated. The possible values are: "mean", "min", "biggest". "biggest" returns the p-value of the biggest connected component.  }
  
  \item{useInteractionSigns}{Logical, should types of interaction be included in the analysis?}
  \item{norm.method}{
Character, type of normalization method to be applied on RNA-Seq data. One of \code{"TMM"} , \code{"DESeq2"} or \code{"none"}. If \code{"none"}, normalization is not performed.  TMM is coupled with voom transformation and DESeq2 with variance stabilizing transformation.
}
}
\section{Warning}{The method is applied witout interaction signs. Argument \code{useInteractionSigns} in the original function \code{testOneGraph} from \code{DEGraph} package is set to \code{FALSE}. }
\value{
A list:
\item{res}{Results from analysis of individual pathways. The first column refers to the overall p-value for a pathway. Then groups of four columns follows. One group refers to one connected component and contains a pair of p-values (without and with Fourier transformation), graph and number of Fourier componets used in the test. The number of groups is equal to the highest number of components in analysed pathways. Components are sorted in the decreasing order of their nodes number.}
 \item{topo.sig }{\code{NULL}, present for the compatibility with outputs from other methods}
 \item{degtest}{A numeric vector of gene-level statistics of all genes in the dataset}
}
\references{
L. Jacob, P. Neuvial, and S. Dudoit. Gains in power from structured two-sample tests of means on graphs. Technical Report arXiv:q-bio/1009.5173v1, arXiv, 2010. 
}
\author{
Ivana Ihnatova
}
\examples{

if (require(DEGraph)) {
  data("Loi2008_DEGraphVignette")

  DEGraph(exprLoi2008, classLoi2008, biocarta[1:10], type="MA")
}
\dontrun{
if (require(gageData)) {

 data(hnrnp.cnts)
 group<-c(rep("sample",4), rep("control",4))
 DEGraph(hnrnp.cnts, group, biocarta[1:10], type="RNASeq", norm.method="TMM")
}
}
}
\keyword{ htest }

