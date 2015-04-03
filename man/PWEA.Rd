\name{PWEA}
\alias{PWEA}

\title{
Function to use PWEA method on microarray or RNA-Seq data
}
\description{
The function runs PWEA method (please see References for the details) on gene expression data matrix, vector specifing to which group a sample belongs and a list of pathway graphs. Briefly, it is a weighted GSEA-like method. The weightes are based on the distance and Pearson's correlation between genes in a pathway.
}
\usage{
PWEA(x, group, pathways, type, preparePaths=TRUE, norm.method=NULL, test.method=NULL, tif=NULL, alpha=0.05, nperm=1000,  ncores=NULL, 
  both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="entrez", convertBy=NULL)
}

\arguments{
  \item{x}{
An \code{ExpressionSet} object or a gene expression data matrix or count matrix, rows refer to genes, columns to samples. Or a list of two data.frames: observed and random (after group permutations) of statistics of differential expression of genes
}
  \item{group}{
Name or number of the phenoData column or a character vector or factor that contains required class assigments
}
  \item{pathways}{
A list of pathways in a form from \code{graphite} package or created by \code{preparePathways()}
}                                                                
  \item{type}{
Type of the data, \code{"MA"} for microarray, \code{"RNASeq"} for RNA-Seq or \code{"DEtable"} for a list of observed and random gene-level statistics
}
 \item{preparePaths}{
Logical, by default the pathways are transformed with \code{preparePathways()}. Use \code{FALSE}, if you have done this transformation separately}
\item{norm.method}{Character, the method to normalize RNAseq data. If \code{NULL} then TMM-normalization is performed. Possible values are: \code{"TMM", "DESeq2", "rLog", "none"}}
\item{test.method}{Character, the method for differentiall expression analysis of RNAseq data. If \code{NULL} then \code{"voomlimma"} is used. Possible values are: \code{"DESeq2", "voomlimma", "vstlimma", "edgeR"}.   }
\item{tif}{A list of Topology Influence Factor's. One slot refers to one pathway. Use \code{prepareTIF()} to create it. It is required only if \code{type=="DEtable"} }
\item{alpha}{Numeric, a theshold value used during TIF calculation}
\item{nperm}{Numeric, number of permutations. Used only if \code{x \%in\% c("MA", "RNASeq")}}
\item{ncores}{Numeric, number of cores. Used only if \code{x \%in\% c("MA", "RNASeq")}. The permutations are calculated in parallel way}
\item{both.directions, maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy}{Arguments for the \code{preparePathways()}}
}

\value{
A list
  \item{res }{A data frame, rows refer to pathways. It contains: Enrichment score for a pathway, p-value and p-value adjusted for multiple hypothesis testing by Benjamini-Hochberg's FDR method. \code{NA}'s if less than 2 nodes are present in the data}
  \item{topo.sig }{A list, topology influence factors for the genes in individual pathways. \code{NULL} if less than 2 nodes are present in the data}
  \item{degtest }{A named vector of statistics from testing the differential expression}

}

\references{

Hung, JH., Whitfield, T. W., Yang, TH., Hu, Z., Weng, Z., DeLisi, Ch. (2010) Identification of functional modules that correlate with phenotypic difference: the influence of network topology, Genome Biology, 11:R23
}
\author{
Ivana Ihnatova
}
\seealso{
 \code{\link{preparePathways}}, \code{\link{prepareTIF}} }
\examples{
\dontrun{
if (require(DEGraph)) {
  data("Loi2008_DEGraphVignette")
  PWEA(exprLoi2008, classLoi2008, biocarta[1:10],  type="MA", nperm=100)
}


if (require(gageData)) {
 data(hnrnp.cnts)
 group<-c(rep("sample",4), rep("control",4))
 PWEA(hnrnp.cnts, group, biocarta[1:10],  type="RNASeq", test="vstlimma", nperm=100)
}
}
}

\keyword{ htest }

