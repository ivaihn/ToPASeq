\name{SPIA}
\alias{SPIA}

\title{
Function to use SPIA method on microarray or RNA-Seq data
}
\description{
The function runs SPIA method on microarray or RNA-Seq data. The implementatio includes the identification of differentially expressed genes and transformation of pathways' topologies to an appropriate form. The SPIA method combines two independent p-values. One p-value comes from overrepresentation analysis and the other is so called pertubation factor. 
}
\usage{
SPIA(x, group, pathways, type="MA", convert=TRUE, IDs = "entrez", 
     gene.stat="logFC", both.directions=TRUE, logFC.th = 2, p.val.th = 0.05, test = NULL, edgeAttrs=defaultEdgeAttrs, ...)
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
  \item{logFC.th}{
A numeric, threshold to select differentially expressed genes according to log Fold-Change. Use negative, if this criterion should not be used.
}
  \item{p.val.th}{
A numeric, threshold to select differentially expressed genes according to p-value from moderated t-test. Use \code{1} if this criterion should not be applied
}
  \item{test}{
Character, which test should be used for differential expression analysis of RNA-Seq data. A character, which approach should be used to calculate differentially expressed genes from RNASeq data. One of \code{"vstlimma"}, \code{"voomlimma"}, \code{"limmaMA"} or \code{"DESeq2"}. \code{"limmaMA"} refers to to \code{limma}-method for microarray data. \code{"vstlimma"} referes to \code{limma} with variance stabilization normalization, \code{"voomlimma"} refers to \code{limma} with \code{voom} transformation. 
}
  \item{edgeAttrs}{Auxiliary data needed for SPIA method. See \code{?edgeAttrs} for more details}
  \item{\dots}{
Preserved for other arguments to \code{spia}. These should not be used.
}
}

\value{
A list:
\item{res}{
A matrix with columns as descibed below:
 pSize - Pathway size, number of genes,
 NDE - Number of differentially expressed genes,
 pNDE - P-value of the overrepresentation part of the method,
 tA - The observed total preturbation accumulation in the pathway,
 pPERT - P-value of the pertubation part of the method,
 p - Combined p-value (overrepresentation and pertubation),
 pFdr - False discovery rate adjusted \code{p},
 pFWER - FWER adjusted \code{p}, 
 Status - If a pathway was identified as Acivated or Inhibited
}
\item{topo.sig}{A list of accumulated pertubation factors and log fold-changes for genes in individual pathways}
\item{degtest}{A numeric vector of gene-level differential expression statistics of all genes in the dataset}
}
\references{
Tarca AL, Draghici S, Khatri P, Hassan SS, Mittal P, Kim JS, Kim CJ, Kusanovic JP, Romero R. A novel signaling pathway impact analysis. Bioinformatics. 2009 Jan 1;25(1):75-82.

Adi L. Tarca, Sorin Draghici, Purvesh Khatri, et. al, A Signaling Pathway Impact Analysis for Microarray Experiments, 2008, Bioinformatics, 2009, 25(1):75-82.

Draghici, S., Khatri, P., Tarca, A.L., Amin, K., Done, A., Voichita, C., Georgescu, C., Romero, R.: A systems biology approach for pathway level analysis. Genome Research, 17, 2007. 
}
\author{
Ivana Ihnatova
}
\examples{

if (require(DEGraph)) {
  data("Loi2008_DEGraphVignette")
  p<-biocarta[1:2]
  SPIA(exprLoi2008, classLoi2008,p, type="MA", logFC.th=-1, IDs="entrez")
}
\dontrun{
if (require(gageData)) {

 data(hnrnp.cnts)
 group<-c(rep("sample",4), rep("control",4))
 SPIA( hnrnp.cnts, group, biocarta[1:10], type="RNASeq",  logFC.th=-1, IDs="entrez", test="vstlimma")
 }
}
}
\keyword{ htest }

