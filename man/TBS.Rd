\name{TBS}
\alias{TBS}

\title{
Function to use TBS method on microarray or RNA-Seq data
}
\description{
A function runs TBS method on a gene expression data matrix or count matrix and vector dividing samples into two groups and a set of pathways from \code{graphite} package. The TBS methods (please see Reference for the details) was adapted to \code{graphite}'s graphs where each node is represented only by one gene.
}
\usage{
TBS(x, group, pathways, type = "MA", test = NULL, convert = TRUE, IDs = "entrez", 
    gene.stat="logFC", both.directions=TRUE, logFC.th = 2, p.val.th = 0.05, nperm = 1000)
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
Type of the data, \code{"MA"} for microarray and \code{"RNA-Seq"} for RNA-Seq
}
  \item{test}{
A character, which approach should be used to calculate differentially expressed genes from RNASeq data. One of \code{"vstlimma"}, \code{"voomlimma"}, \code{"limmaMA"} or \code{"DESeq2"}. \code{"limmaMA"} refers to to \code{limma}-method for microarray data. \code{"vstlimma"} referes to \code{limma} with variance stabilization normalization, \code{"voomlimma"} refers to \code{limma} with \code{voom} transformation
}
 \item{convert}{
Logical, should nodes' labels be converted?
}
  \item{IDs}{
Character, which gene IDs should be used. One of \code{"entrez"} or \code{"symbol"}
}
  \item{gene.stat}{Character, which gene-level statistic should be exported? \code{"logFC"} stands for log fold-change and \code{"stats"} for the modified t-statistic from \code{limma} package}
   \item{both.directions}{Logical, should undirected edges be substitued by two directed edges (one in each direction) or only directed edges are considered NOTE: this is different from the other methods implemented in this package}
  \item{logFC.th}{
Numeric, threshold for log fold-change of a gene to identify the gene as differentially expressed. Use negative if you don't want any threshold to be applied
}
  \item{p.val.th}{
Numeric, threshold for p-values of moderated t-test (from \code{limma} package) of a gene to identify the gene as significantly differentially expressed. Use \code{1} if you don't want any threshold to be applied
}
  \item{nperm}{
Numeric, number of permutations
}
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
\examples{

if (require(DEGraph)) {
  data("Loi2008_DEGraphVignette")
  TBS( exprLoi2008, classLoi2008, biocarta[1:10], type="MA",  logFC.th=-1, nperm=100)
}
\dontrun{
if (require(gageData)) {

 data(hnrnp.cnts)
 group<-c(rep("sample",4), rep("control",4))
 
 TBS(hnrnp.cnts, group, biocarta[1:10], type="RNASeq", logFC.th=-1, nperm=100, test="vstlimma")
}
}
}

\keyword{ htest }

