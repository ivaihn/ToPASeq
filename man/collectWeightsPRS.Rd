\name{collectWeightsPRS}
\alias{collectWeightsPRS}
\alias{collectWeightsSPIA}
\alias{prepareTIF}

\title{
Function to calculate gene-level weights for topology-based pathway analysis
}
\description{
The functions calculate gene-level weigths defined in various topology-based pathway analysis methods (PRS, SPIA, PWEA). In PRS, it is the number of downstream differentially expressed genes. TIF, the statistic defined in PWEA, is related to the ratio of correlation and distance of genes. SPIA defines the so called net pertubation factors. 
}
\usage{
collectWeightsPRS(de, all, pathways)
collectWeightsSPIA(de, all, pathways)
prepareTIF(pathways, exprs, alpha)
}

\arguments{
  \item{de}{
Named numeric vector, the log fold-changes of the differentially expressed genes
}
  \item{all}{
Character vector of all genes meassured in the experiment
}
  \item{pathways}{
A list of pathways, each pathway is an object of class \code{Pathway} transformed via \code{preparePathways()} for the particular method}

\item{exprs}{A numeric matrix, gene expression data matrix, rows refer to genes, columns to samples}
\item{alpha}{Numeric, a threshold to control the magnitude. In TIF calculation, the effect of a gene on a few nearby and tightly correlated genes can be washed out if the gene influences many other genes weakly. The threshold supresses this washing-out
}
}

\value{
A list, each slot is a vector of gene-level weights for one pathway
}

\author{
Ivana Ihnatova  
}
\examples{
pathways<-pathways("hsapiens","kegg")[1:3]
de<-setNames(rnorm(30),sample(nodes(pathways[[1]]),30))
all<-nodes(pathways[[1]])

path<-preparePathways(pathways[1:3], method="SPIA", genes=all, both.direction=TRUE, convertTo="none")
collectWeightsSPIA(de, all, path)
}
\keyword{ univar }

