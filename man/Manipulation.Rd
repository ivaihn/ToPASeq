\name{Manipulation}
\docType{methods}
\alias{Manipulation}
\alias{intersection,pathway,pathway-method}
\alias{join,pathway,pathway-method}
\alias{union,pathway,pathway-method}
\alias{subGraph,character,pathway-method}

\alias{clearNode,character,pathway-method}
\alias{removeEdge,character,character,pathway-method}
\alias{removeNode,character,pathway-method} 
\alias{addEdge,character,character,pathway,numeric-method}
\alias{addNode,character,pathway,list-method}
\alias{nodes<-,pathway,character-method}
\alias{changeDirection,character,character,pathway,character-method} 
\alias{changeDirection}
\alias{changeInteraction,character,character,pathway,character,logical-method}
\alias{changeInteraction}

\title{
Methods for manual editing of \code{pathway} class from \pkg{graphite}
}
\section{Methods}{
\describe{
All of the methods below return an object of class \code{pathway} with modified topology.

\item{intersection \code{signature(x = "pathway", y = "pathway")}}{compute the intersection of the two supplied graphs. They must have identical nodes.}

\item{join \code{signature(x = "pathway", y = "pathway")}}{returns the joining of the two graphs. It is similar to \code{intersection} but does not require the identical nodes}

\item{union \code{signature(x = "pathway", y = "pathway")}}{compute the union of the two supplied graphs. They must have identical nodes.}

\item{subGraph \code{signature(snodes = "character", graph = "pathway")}}{Given a set of nodes and a pathway this function creates and returns subgraph with only the supplied nodes and any edges between them  }          

\item{clearNode \code{signature(node = "character", object = "pathway")}}{Clears all edges incoming and outgoing edges from \code{node}}

\item{removeEdge \code{signature(from = "character", to = "character", graph = "character")}}{removes all directed edges starting in \code{from} and ending in \code{to} and undirected edges between \code{from} and \code{to} }

\item{removeNode \code{signature(node = "character", object = "pathway")}}{removes node(s) \code{node} from a pathway \code{object}}

\item{nodes<- \code{signature(object = "pathway", value = "character")}}{sets node labels of pathway \code{object} to \code{value}}

}
}
\description{ \code{changeInteraction} and \code{changeDirection} are a new generic function designed for \code{pathway} class}

\value{
A \code{pathway}
}
\author{
Ivana Ihnatova
}

\keyword{ methods }

