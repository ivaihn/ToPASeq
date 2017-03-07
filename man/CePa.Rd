\name{CePa}
\alias{CePa}

\title{
Function to use CePa method on microarray or RNA-Seq data
}
\description{
The function runs CePa method on microarray or RNA-Seq data. The implementatio includes the identification of differentially expressed genes and transformation of pathways' topologies to an appropriate form. The ORA version of the CePa method is implemented and covers centralities: equal-weight, in-degree, out-degree, in-reach, out-reach and betweenness.
}
\usage{
CePa(x, group, pathways, type, preparePaths=TRUE, norm.method=NULL, test.method=NULL, p.th=0.05, logFC.th=2, nperm=1000, 
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
A list:
\item{res}{
A matrix, each row refers to one pathway, each column to one centrality and the value is a p-value. 
}
\item{topo.sig}{A list of weights for genes (nodes) in individual pathways}
\item{degtest}{A numeric vector of gene-level differential expression statistics of all genes in the dataset}
}
\references{
Gu Z., Liu J., Cao K., Zhang J., Wang J.: Centrality-based pathway enrichment: a systematic approach for finding significant pathways dominated by key genes. BMC Systems Biology 2012, 6:56 
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
  CePa(exprLoi2008, classLoi2008,pathways, type="MA", logFC.th=-1)
}
\dontrun{
if (require(gageData)) {

 data(hnrnp.cnts)
  hnrnp.cnts<-hnrnp.cnts[rowSums(hnrnp.cnts)>0,]
 group<-c(rep("sample",4), rep("control",4))
 
  pathways<-pathways("hsapiens","biocarta")[1:10]
 SPIA( hnrnp.cnts, group, pathways, type="RNASeq",  logFC.th=-1, IDs="entrez", test="vstlimma")
 }
}
}
\keyword{ htest }

