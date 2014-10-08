\name{TopologyGSA}
\alias{TopologyGSA}

\title{
Function to use TopologyGSA method on microarray or RNA-Seq data
}
\description{
TopologyGSA method uses graphical models to test the differential expression of a pathway. It also highlights pathway componenets involved in the deregulation. 
}
\usage{
TopologyGSA(x, group, pathways, type = "MA", convert = TRUE, IDs = "entrez",  both.directions=TRUE,
            test="mean", testCliques=FALSE, alpha=0.05, nperm = 10000, 
            norm.method = NULL, maxNodes = 150)
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
 \item{both.directions}{Logical, should undirected edges be substitued by two directed edges (one in each direction) or only by one directed edge oriented as in \code{edges()}}
  \item{test}{
Either \code{"var"} and \code{"mean"}. Determine the type of test used by topologyGSA.
}
\item{testCliques}{Logical, should also be individual cliques for each pathway tested? Please, be aware that this calculation may take quite a lot of time}
  \item{alpha}{
Numeric, threshold for statistical significance
}
  \item{nperm}{
Numeric, Number of permutations
}
  \item{norm.method}{
Character, type of normalization method to be applied on RNA-Seq data. One of \code{"TMM"} , \code{"DESeq2"} or \code{"none"}. If \code{"none"}, normalization is not performed. TMM is coupled with voom transformation and DESeq2 with variance stabilizing transformation.
}
\item{maxNodes}{Pathways with more nodes than specified value will be omitted from the analysis}

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
  p<-biocarta[1:10]
  TopologyGSA(exprLoi2008, classLoi2008, p, type="MA", test="mean", alpha=0.05, nperm=200)
  }
  
if (require(gageData)) {

 data(hnrnp.cnts)
 group<-c(rep("sample",4), rep("control",4))
 TopologyGSA(hnrnp.cnts, group, biocarta[1:10], type="RNASeq",test="mean", alpha=0.05, 
   nperm=200, norm.method="TMM")
}
}
}
\keyword{ htest }

