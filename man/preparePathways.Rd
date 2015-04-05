\name{preparePathways}
\alias{preparePathways}

\title{
Function to prepare pathways for topology-based pathway analysis
}
\description{
Functions transforms pathways from \code{graphite} package (stored as \code{pathway-class}) into formats required in the particular topology-based method implemented in this package. It also converts identifiers in the pathways and filters pathways according to several criteria.
}
\usage{
preparePathways(pathways, method, both.directions, genes, maxNodes = 150, minEdges = 0, commonTh = 2, filterSPIA = FALSE, convertTo = "entrez", convertBy = NULL, EdgeAttrs = NULL)
}

\arguments{
  \item{pathways}{
A list of pathways, individual pathways are objcets of class \code{pathway}
}
  \item{method}{
A character, the pathways will be transformed according to the needs of the particular method. Possible values are: \code{"TAPPA", "PRS", "PWEA", "TopologyGSA", "clipper", "DEGraph","SPIA" }
}
  \item{both.directions}{
Logical, indicates how should be the undirected edges directed. If \code{TRUE}, an undirected edge is substitued with two directed edges with opposite directions (e.g. A-B becomes A->B and B->A). If \code{FALSE}, then an undirected edge is substitued with one directed edge which preserves the order of nodes (e.g. A-B becomes A->B).
}
  \item{genes}{
Character vector, vector of gene identifiers in the expression data
}
  \item{maxNodes}{
Numeric, maximal number of nodes. Pathways with more nodes are filtered out.
}
  \item{minEdges}{
Numeric, minimal number of edges. Pathways with less edges are filtered out.
}
  \item{commonTh}{
Numeric, threshold for number of nodes present in the data. Pathways with less node-identifiers matching to \code{genes} are filtered out.
}
  \item{filterSPIA}{
Logical, if \code{TRUE} applies filter defined in the SPIA method (relates to the calculation of inversion matrix).
}
  \item{convertTo}{
Character, one of the \code{"entrez"}, \code{"symbol"} or \code{"none"}. If \code{"none"} no conversion is performed. Otherwise, the function converts node-identifiers in pathways from \code{graphite} into either EntrezIDs or gene symbols. Uses \code{org.Hs.eg.db} package for the mapping. 
}
  \item{convertBy}{
Named character vector, names of the elemenet must match the node-identifiers and the values are the new identifiers to be replaced. This is a more general option designed for pathways outside \code{graphite}.  
}
  \item{EdgeAttrs}{
A list of two tables required for the filter from SPIA method. See \code{makeDefaultEdgeData} for the details.
}
}

\value{
A list of the transformed pathways}

\author{
Ivana Ihnatova
}

\seealso{
 \code{\link{makeDefaultEdgeData}}}
\examples{
#Creating dummy set of genes
set.seed(123)
genes<-unname(unlist(lapply(kegg[1:3], nodes)))
genes<-sample(genes, length(genes)*0.9)

#Applying the function
paths<-preparePathways(kegg[1:3], "TAPPA", TRUE, genes, maxNodes=65, convertTo="none")
paths


}
\keyword{manip}

