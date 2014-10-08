\name{TAPPA}
\alias{TAPPA}

\title{
Function to use TAPPA method on microarray or RNA-Seq data
}
\description{
The functions analyses the differential expression of pathways via TAPPA method. Expression is compared between two groups of samples by Mann-Whitney test. P-values are later adjusted for multiple hypothesis testing by Benjamini-Hochberg's FDR method.
}
\usage{
TAPPA(x, group, pathways, type = "MA", convert = TRUE, IDs = "entrez", 
      gene.stat="logFC", both.directions=TRUE, normalize = TRUE, verbose = FALSE, norm.method = NULL)
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
  \item{gene.stat}{Character, which gene-level statistic should be exported? \code{"logFC"} stands for log fold-change and \code{"stats"} for the modified t-statistic from \code{limma} package}  ]
 \item{both.directions}{Logical, should undirected edges be substitued by two directed edges (one in each direction) or only by one directed edge oriented as in \code{edges()}}
  \item{normalize}{
Logical, should data be normalized? 
}
  \item{verbose}{
Logical, if \code{TRUE} names of the pathways are printed as they are analysed
}
  \item{norm.method}{
Character, type of normalization method to be applied on RNA-Seq data. One of \code{"TMM"} , \code{"DESeq2"} or \code{"none"}. If \code{"none"}, normalization is not performed. TMM is coupled with voom transformation and DESeq2 with variance stabilizing transformation.
}
}

\value{ 
A list,
\item{res}{A data frame, rows refer to pathways. Columns contain: number of valid PCI-scores, median, min and max of the PCI scores for each group of samples (denoted with \code{1} or \code{2}), p-value of the Mann-Whitney test (\code{p.val}) and adjusted p-value (\code{p.adj}). If less than two nodes are present in the data, the function puts \code{NA}'s in all columns. }
\item{topo.sig}{\code{NULL}, it is preserved for the compatibility with other methods implemented in this package}
\item{degtest}{A numeric vector of gene-level differential expression statistics}
}
\references{
Gao, S. and Wang, X. (2007) TAPPA: topological analysis of pathway phenotype association. Bioinformatics, 23, pages 3100-3102
}
\author{Ivana Ihnatova}

\examples{

if (require(DEGraph)) {
  data("Loi2008_DEGraphVignette")
  TAPPA(exprLoi2008, classLoi2008, biocarta[1:10], type="MA")
}

\dontrun{
if (require(gageData)) {

 data(hnrnp.cnts)
 group<-c(rep("sample",4), rep("control",4))
 TAPPA( hnrnp.cnts, group, biocarta[1:10], type="RNASeq", norm.method="TMM")
}
}
}


\keyword{ htest }

