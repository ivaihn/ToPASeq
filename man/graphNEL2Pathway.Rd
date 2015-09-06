\name{graphNEL2Pathway}
\alias{graphNEL2Pathway}

\title{
Function to coerce a graphNEL to a Pathway
}
\description{
The function coerces a graphNEL to a Pathway. It attempts to recover the edge types from \code{"edgeType"} attribute of \code{edgeData}. The result contains only the edge types present in the graph. If the \code{edgeData} do not contain this attribute, then \code{"process(indirect effect)"} is used in order to preserve directionality. 
}
\usage{
graphNEL2Pathway(graph, name = "pathway", ident = "unknown", database = "unknown", species = "unknown", date = NULL)
}

\arguments{
  \item{graph}{
A \code{graphNEL} object to be coerced.
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
A coerced \code{Pathway}
}
\note{
When this function is applied on \code{x} as reversed operation to \code{pathwayGraph} then the order of the edges may differ as well as the directionality of \code{"process(indirect)"} edges as they are set as \code{undirected} by graphNEL2Pathway.
}
\author{
Ivana Ihnatova
}


\examples{
pathway<-pathways("hsapiens","kegg")[[1]]
pathway<-pathwayGraph(pathway)
pathway
graphNEL2Pathway(pathway)

set.seed(123)
rg <- randomEGraph(LETTERS[1:20], edges = 30)
p<-graphNEL2Pathway(rg)
p
head(edges(p))
}

\keyword{manip}

