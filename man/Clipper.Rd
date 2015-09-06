\name{clipper}
\alias{clipper}

\title{
Function to use clipper method on microarray or RNA-Seq data
}
\description{
clipper is a method for topological gene set analysis. It implements a two-step empirical approach based on the exploitation of graph decomposition into a junction tree to reconstruct the most relevant signal path. In the first step clipper selects significant pathways according to statistical tests on the means and the concentration matrices of the graphs derived from pathway topologies. Then, it "clips" the whole pathway identifying the signal paths having the greatest association with a specific phenotype. 
}
\usage{
clipper(x, group, pathways, type, preparePaths=TRUE, norm.method=NULL, test.method=NULL,  method="mean", testCliques=FALSE, nperm=1000, alphaV=0.05, b=1000, permute=TRUE,
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
Type of the input data, \code{"MA"} for microarray and \code{"RNASeq"} for RNA-Seq
}
\item{preparePaths}{
Logical, by default the pathways are transformed with \code{preparePathways()}. Use \code{FALSE}, if you have done this transformation separately}
\item{norm.method}{Character, the method to normalize RNAseq data. If \code{NULL} then TMM-normalization is performed. Possible values are: \code{"TMM", "DESeq2", "rLog", "none"}}
\item{test.method}{Character, the method for differentiall expression analysis of RNAseq data. If \code{NULL} then \code{"voomlimma"} is used. Possible values are: \code{"DESeq2", "voomlimma", "vstlimma", "edgeR"}. This analysis is needed only for the visualization.   }
\item{method}{Character, \code{"mean"}  or \code{"var"}, the kind of test to perform on the cliques}
\item{testCliques}{Logical, if \code{TRUE} then the test is applied also on the cliques of the each pathway. It is a very time consuming calculation, especially for many or big pathways}
\item{nperm}{Number of permutations}
\item{alphaV}{Numeric, the threshold for variance test. The calculation of mean test depends on the result of variance test. }
\item{b}{number of permutations for mean analysis}
\item{permute}{always performs permutations in the concentration matrix test. If FALSE, the test is made using the asymptotic distribution of the log-likelihood ratio. This option should be use only if samples size is >=40 per class}
\item{both.directions, maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy}{Arguments for the \code{preparePathways()}}
 
}
\value{
A list,
\item{res}{A list. First slot is a data frame containing p-values and q-values of mean and variance tests on pathways. The second slot is a list containing data.frames of the most affected paths in each pathway. The columns of the data frames contain: 1 - Index of the starting clique 2 - Index of the ending clique 3 - Index of the clique where the maximum value is reached 4 - length of the path 5 - maximum score of the path 6 - average score along the path 7 - percentage of path activation 8 - impact of the path on the entire pathway 9 - clique involved and significant 10 - clique forming the path 11 - genes forming the significant cliques 12 - genes forming the path}
\item{topo.sig}{if \code{testCliques=TRUE}, a list where each slot contains the pvalues and a list of cliques in one pathway. \code{NULL} otherwise}
\item{degtest}{A data.frame of gene-level differential expression statistics}
}
\references{
Martini P, Sales G, Massa MS, Chiogna M, Romualdi C. Along signal paths: an empirical gene set approach exploiting pathway topology. Nucleic Acids Res. 2013 Jan 7;41(1):e19. doi: 10.1093/nar/gks866. Epub 2012 Sep 21. PubMed PMID: 23002139; PubMed Central PMCID: PMC3592432. 
}
\author{
Ivana Ihnatova
}
\note{
If there are NA's only in columns 3 to 7, then a junction tree could not be formed. 
}
\seealso{
 \code{\link{preparePathways}}}
\examples{


if (require(DEGraph)) {
  data("Loi2008_DEGraphVignette")
 pathways<-pathways("hsapiens","kegg")[1]
  clipper( exprLoi2008, classLoi2008, pathways,type="MA", convertTo="none")
}   

\dontrun{
if (require(gageData)) {

 data(hnrnp.cnts)
 hnrnp.cnts<-hnrnp.cnts[rowSums(hnrnp.cnts)>0,]
 group<-c(rep("sample",4), rep("control",4))
  pathways<-pathways("hsapiens","kegg")[1:3]
 clipper(hnrnp.cnts, group,pathways, type="RNASeq",  norm.method="TMM", convertTo="none")
 }
}
    
}
\keyword{htest }

