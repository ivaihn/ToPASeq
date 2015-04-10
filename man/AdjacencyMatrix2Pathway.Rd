\name{AdjacencyMatrix2pathway}
\alias{AdjacencyMatrix2pathway}

\title{
Function to coerce an adjacency matrix to a pathway
}
\description{
The function coerces an adjacency matrix to a pathway. Two types of matrices are allowed. The first one, where \code{1} denotes an edge between two nodes and \code{0} otherwise. This matrix is coerced into a simply pathway were type of all edges is set to \code{"process"}. The second type of adjacency matrix contains: \code{1} for an activation, \code{-1} for an inhibition and \code{0} otherwise (=no edge between two nodes). In this case, activations are set to \code{"process(activation)"} and inhibition to \code{"process(inhibition)"}. The symetricity of the matrix is used to decide between directed and undirected graph. Symmetric matrix is expected for undirected graph and only the lower triangle of the matrix is used to extract the edges of the graph.
}

\usage{
AdjacencyMatrix2pathway(adjmat, name = "pathway", ident = "unknown", database = "unknown", species = "unknown", date = NULL)
}
\arguments{
  \item{adjmat}{
An adjacency matrix describing the pathway topology
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
An object of class \code{pathway}}
\author{
Ivana Ihnatova
}


\examples{
genes<-paste("gene", 1:10, sep="")
adjmat<-matrix(sample(c(0,0,0,0,1), 100, TRUE),10,10, dimnames=list(genes,genes))
p<-AdjacencyMatrix2pathway(adjmat)
head(edges(p))

adjmat<-matrix(sample(c(0,0,0,0,1,-1), 100, TRUE),10,10, dimnames=list(genes,genes))
p<-AdjacencyMatrix2pathway(adjmat)
head(edges(p))

}

\keyword{ manip }

