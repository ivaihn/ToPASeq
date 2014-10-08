\name{PWEA}
\alias{PWEA}

\title{
Function to use PWEA method on microarray or RNA-Seq data
}
\description{
The function runs PWEA method (please see References for the details) on gene expression data matrix, vector specifing to which group a sample belongs and a list of pathway graphs. Briefly, it is a weighted GSEA-like method. The weightes are based on the distance and Pearson's correlation between genes in a pathway.
}
\usage{
PWEA(x, group, pathways, type = "MA", test=NULL, convert = TRUE, IDs = "entrez",
     gene.stat="logFC", both.directions=TRUE, alpha = 0.05, nperm = 5000)
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
  \item{test}{
A character, which approach should be used to calculate differentially expressed genes from RNASeq data. One of \code{"vstlimma"}, \code{"voomlimma"}, \code{"limmaMA"} or \code{"DESeq2"}. \code{"limmaMA"} refers to to \code{limma}-method for microarray data. \code{"vstlimma"} referes to \code{limma} with variance stabilization normalization, \code{"voomlimma"} refers to \code{limma} with \code{voom} transformation, \code{DESeq2} method takes much longer and therefore it is not recommended.}
  \item{convert}{
Logical, should nodes' labels be converted?
}
  \item{IDs}{
Character, which gene IDs should be used. One of \code{"entrez"} or \code{"symbol"}
}
  \item{gene.stat}{Character, which gene-level statistic should be exported? \code{"logFC"} stands for log fold-change and \code{"stats"} for the modified t-statistic from \code{limma} package   }
   \item{both.directions}{Logical, should undirected edges be substitued by two directed edges (one in each direction) or only by one directed edge oriented as in \code{edges()}}
  \item{alpha}{
Numeric, a parameter used when calculating genes weights
}
  \item{nperm}{
Numeric, number of permutations
}
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

