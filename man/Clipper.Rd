\name{Clipper}
\alias{Clipper}

\title{
Function to use clipper method on microarray or RNA-Seq data
}
\description{
clipper is a method for topological gene set analysis. It implements a two-step empirical approach based on the exploitation of graph decomposition into a junction tree to reconstruct the most relevant signal path. In the first step clipper selects significant pathways according to statistical tests on the means and the concentration matrices of the graphs derived from pathway topologies. Then, it "clips" the whole pathway identifying the signal paths having the greatest association with a specific phenotype. 
}
\usage{
Clipper(x, group, pathways, type = "MA", convert = TRUE, IDs = "entrez",  both.directions=TRUE,
        test="mean", testCliques=FALSE, nperm=1000, norm.method = NULL)
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
Character, a test to perform on the cliques. It can be mean or variance
}
\item{testCliques}{Logical, should also be individual cliques for each pathway tested? Please, be aware that this calculation may take quite a lot of time}
\item{nperm}{Numeric, number of permutations in tests}
  \item{norm.method}{
Character, type of normalization method to be applied on RNA-Seq data. One of \code{"TMM"}, \code{"DESeq2"} or \code{"none"}. If \code{"none"}, normalization is not performed.    TMM is coupled with voom transformation and DESeq2 with variance stabilizing transformation.
}
}
\value{
A list,
\item{res}{A data frame, where rows refer to pathways and columns are organized as follows: 1 - P-value from comparing the covariance matricesof whole pathway 2 - P-value from comparing the mean expression of whole pathway 3 - maximum score of the paths 4 - percentage of activation of the maximal scoring path 5 - impact of the maximal scoring path on the entire pathway 6 - genes forming the significant cliques 7 - genes forming the maximal scoring path)}
\item{topo.sig}{if \code{testCliques=TRUE}, a list where each slot contains the pvalues and a list of cliques in one pathway. \code{NULL} otherwise}
\item{degtest}{A numeric vector of gene-level differential expression statistics}
}
\references{
Martini P, Sales G, Massa MS, Chiogna M, Romualdi C. Along signal paths: an empirical gene set approach exploiting pathway topology. Nucleic Acids Res. 2013 Jan 7;41(1):e19. doi: 10.1093/nar/gks866. Epub 2012 Sep 21. PubMed PMID: 23002139; PubMed Central PMCID: PMC3592432. 
}
\author{
Ivana Ihnatova
}
\note{
If only NA's are returned for a pathway, then less than 2 genes are present in the data. If there are NA's only in columns 3 to 7, then a junction tree could not be formed. 
}
\examples{

if (require(clipper)) {
if (require(DEGraph)) {
  data("Loi2008_DEGraphVignette")
  Clipper( exprLoi2008, classLoi2008, kegg[1],type="MA")
}
}
\dontrun{
if (require(gageData)) {

 data(hnrnp.cnts)
 group<-c(rep("sample",4), rep("control",4))
 Clipper(hnrnp.cnts, group,kegg[1:10], type="RNASeq",  norm.method="TMM")
 }
 }
    
}
\keyword{htest }

