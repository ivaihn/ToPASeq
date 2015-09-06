\name{TopologyGSA}
\alias{TopologyGSA}

\title{
Function to use TopologyGSA method on microarray or RNA-Seq data
}
\description{
TopologyGSA method uses graphical models to test the differential expression of a pathway. It also highlights pathway componenets involved in the deregulation. 
}
\usage{
TopologyGSA(x, group, pathways, type, preparePaths=TRUE, norm.method=NULL, test.method=NULL , method="mean", alpha=0.05, testCliques=FALSE, ..., 
   both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL )
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
  \item{method}{
Either \code{"var"} and \code{"mean"}. Determine the type of test used by topologyGSA.
}
 \item{alpha}{
Numeric, threshold for statistical significance of variance test. It influences the method for the mean test
}
\item{testCliques}{Logical, if \code{TRUE}, then the test is also performed on individual cliques. It can be very computationally complex.}
\item{...}{Other arguments to be passed to the method. See details for better explanation}
\item{both.directions, maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy}{Arguments for the \code{preparePathways()}}
 
}

\details{
The method requires a Directed Acyclic Graph (DAG). Therefore if a pathway contain also undirected or bidirected edges and error is thrown. 

The user can further specify for the mean test:
\enumerate{
\item \strong{perms} number of permutations of the test,
\item \strong{paired}logical flag. If TRUE Hotelling test for paired samples is calculated and the test on the variances is not performed 
}

Or for the variance test:
\enumerate{
\item \strong{variance}logical flag. If TRUE the estimates of the covariance matrices are included in the result.
\item \strong{s1}First group covariance matrix estimation.
\item \strong{s2}Second group covariance matrix estimation.
}
}
\value{
A list
\item{res}{a list with one entry for each successfully analyzed pathway }

\item{topo.sig}{if \code{testCliques=TRUE}, a list where each slot contains the pvalues and a list of cliques in one pathway. \code{NULL} otherwise}
\item{degtest}{A numeric vector of gene-level differential expression statistics}

}
\references{
Massa MS, Chiogna M, Romualdi C. Gene set analysis exploiting the topology of a pathway. BMC System Biol. 2010 Sep 1;4:121.
}
\author{
Ivana Ihnatova
}

\examples{
\dontrun{
if (require(DEGraph)) {
  data("Loi2008_DEGraphVignette")
  pathways<-pathways("hsapiens","biocarta")[1:10]

  
  TopologyGSA(exprLoi2008, classLoi2008, pathways, type="MA", method="mean", alpha=0.05, perms=200)
    TopologyGSA(exprLoi2008, classLoi2008, pathways, type="MA", method="mean", alpha=0.05, perms=200, testCliques=TRUE)
  }
  
if (require(gageData)) {

 data(hnrnp.cnts)
 group<-c(rep("sample",4), rep("control",4))
 hnrnp.cnts<-hnrnp.cnts[rowSums(hnrnp.cnts)>0,]
 pathways<-pathways("hsapiens","biocarta")[1:10]
 TopologyGSA(hnrnp.cnts, group,pathways, type="RNASeq",method="mean", alpha=0.05, 
   perms=200, norm.method="TMM")
}
}
}
\keyword{ htest }

