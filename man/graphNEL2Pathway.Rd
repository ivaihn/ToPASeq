\name{graphNEL2pathway}
\alias{graphNEL2pathway}

\title{
Function to coerce a graphNEL to a pathway
}
\description{
The function coerces a graphNEL to a pathway. It attempts to recover the edge types from \code{"edgeType"} attribute of \code{edgeData}. The result contains only the edge types present in the graph. If the \code{edgeData} do not contain this attribute, then \code{"process(indirect effect)"} is used in order to preserve directionality. 
}
\usage{
graphNEL2pathway(graph, name = "pathway", ident = "unknown", database = "unknown", species = "unknown", date = NULL)
}

\arguments{
  \item{graph}{
A graphNEL object to be coerced.
}
  \item{name}{
A character, name of the pathway. Defaults to \code{"pathway"}
}
  \item{ident}{
A character, type of the identificators, e.g \code{"gene symbol"}
}
  \item{database}{
A character, the name of the database the topology comes from
}
  \item{species}{
A character, the species to which the topology belong
}
  \item{date}{
A date, the date the topology was created
}
}


\value{
A coerced \code{pathway}
}
\author{
Ivana Ihnatova
}



\seealso{
 \code{\link[graphite]{pathway-class}}, \code{\link[graphite]{pathwayGraph}}
}
\examples{
kegg[[1]]
pathwayGraph(kegg[[1]])
graphNEL2pathway(pathwayGraph(kegg[[1]]))

set.seed(123)
rg <- randomEGraph(LETTERS[1:20], edges = 30)
p<-graphNEL2pathway(rg)
p
head(edges(p))
}

\keyword{manip}

