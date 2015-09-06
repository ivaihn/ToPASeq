\name{convertIdentifiersByVector}
\alias{convertIdentifiersByVector}

\title{
Function to convert identifiers in pathways by user specified vector
}
\description{
The function converts identifiers of nodes in a pathway. It uses the user specified named vector for the conversion. 
}
\usage{
convertIdentifiersByVector(pathway, conv.table, id.type="unknown")
}

\arguments{
  \item{pathway}{
An object of class \code{Pathway}
}
  \item{conv.table}{
A named vector in which names correspond to the identifiers present in the pathway and values are the new identifiers to which conversion happens
}
  \item{id.type}{
A character, the type of the identifiers provided e.g \code{"TAIR"} for TAIR numbers. This is for informative purposes only.
}
}

\value{
A \code{Pathway} in which identifiers have been converted
}

\author{
Ivana Ihnatova
}

\seealso{
 \code{\link[graphite]{convertIdentifiers}}     }

\examples{
g<-kegg[["Asthma"]]
conv<-setNames(paste("gene", 1:length(nodes(g)), sep=""), nodes(g))
gc<-convertIdentifiersByVector(g, conv, "dummy")
nodes(gc)
edges(gc)
}
\keyword{manip}

