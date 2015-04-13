\name{estimateCF}
\alias{estimateCF}

\title{
Function to estimate multi-subunit protein complexes and gene families in a pathway
}
\description{
Function estimates the multi-subunit protein complexes and gene families in a pathway. A protein complex consists of proteins connected by undirected binding interaction. A gene family is a set of nodes with same outgoing and/or incomig edges. 
}
\usage{
estimateCF(graph)
}

\arguments{
  \item{graph}{
An object of class \code{pathway}
}
}

\value{
  \item{complexes }{A list of estimated protein complexes'}
  \item{famillies }{A list of estimated gene famillies}
The function attempts to assign a representative name to each gene family. The representative name is a common part of the names of individual genes. This approach, however, may lead to ambiguities or missings. Then a general name in a form of \code{family1}, \code{family2}, etc. All the complexes are named analogously as \code{complex1}, \code{complex2}.
}

\author{
Ivana Ihnatova
}


\seealso{
 \code{\link{reduceGraph}}
}
\examples{
pathway<-as(kegg[[1]], "pathway")
estimateCF(pathway)
}
\keyword{manip }

