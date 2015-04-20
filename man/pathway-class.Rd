\name{pathway-class}
\Rdversion{1.1}
\docType{class}
\alias{pathway-class}
\alias{edges,pathway,ANY-method}
\alias{nodes,pathway-method}
\alias{show,pathway-method}

\alias{intersection,pathway,pathway-method}
\alias{join,pathway,pathway-method}
\alias{union,pathway,pathway-method}
\alias{subGraph,character,pathway-method}
\alias{degree,pathway,character-method}
\alias{degree,pathway,missing-method}
\alias{numNoEdges,pathway-method}
\alias{mostEdges,pathway-method}
\alias{acc,pathway,character-method}
\alias{connComp,pathway-method}
\alias{edges,pathway,character-method}
\alias{isAdjacent,pathway,character,character-method}
\alias{isConnected,pathway-method}
\alias{isDirected,pathway-method}
\alias{edgemode,pathway-method}
\alias{numEdges,pathway-method}
\alias{numNodes,pathway-method}
\alias{edgeNames,pathway-method}

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
\alias{convertIdentifiers,pathway,character-method}

\title{Class \code{"pathway"}}
\description{
  This class represents a biological pathway. It is a class from the package \code{graphite} up to version \code{1.12}.  \code{changeInteraction} and \code{changeDirection} are a new generic function designed for \code{pathway} class
}
\section{Slots}{
  \describe{
    \item{\code{title}:}{Object of class \code{"vector"}.}
    \item{\code{nodes}:}{Object of class \code{"vector"}.}
    \item{\code{edges}:}{Object of class \code{"data.frame"}.
      The \code{data.frame} has one row for each edge and four columns:
      source node, destination node, edge direction and type.
    }
    \item{\code{ident}:}{Object of class \code{"vector"}.
      The type of identifiers used to label the nodes. 

    }
    \item{\code{database}:}{Object of class \code{"vector"}.
      The name of the database from which the pathway was retrieved.
    }
    \item{\code{species}:}{Object of class \code{"character"}.
      The name of the species the pathway belongs to.
    }
    \item{\code{timestamp}:}{Object of class \code{"Date"}.
      The date when the pathway was retrieved.
    }
  }
}
\section{Methods}{
  \describe{
    \item{edges}{\code{signature(object = "pathway")}:
      retrieves the \code{data.frame} describing the pathway edges.
    }

    \item{nodes}{\code{signature(object = "pathway")}: retrieves the
      \code{vector} enumerating the identifiers of the pathway nodes.
    }
  }
  The methods below perfom basic topological analysis of a pathway. They were defined as generic in \code{graph} for \code{graph} class. They were implemented for \code{pathway} in this package
  \describe{

\item{degree \code{signature(object = "pathway", Nodes = "character")}}{Returns the number of incoming or outgoing edges for nodes in \code{Nodes} }
\item{degree \code{signature(object = "pathway", Nodes = "missing")}}{Returns the number of incoming or outgoing edges for all nodes in \code{object} }
\item{numNoEdges \code{signature(objGraph = "pathway")}}{Returns the number of nodes without any edge}
\item{mostEdges \code{signature(objGrapg = "pathway")}}{Returns the nodes with most edges}
\item{acc \code{signature(object = "pathway", index = "character")}}{Returns the set of nodes accessible from nodes in \code{index}. The undirected edges are considered as bidirected (directed in both directions)}
\item{connComp \code{signature(object = "pathway")}}{Returns the connected components present in a pathway. They are returned as list where each slot refers to one component and contains the relevant nodes. The undirected edges are considered as bidirected (directed in both directions)}
\item{edges \code{signature(object = "pathway", which = "character")}}{Returns the edges relevant to node(s) in \code{which}}
\item{isAdjacent \code{signature(object = "pathway", from = "character", to = "character")}}{Returns whether nodes in \code{from} and \code{to} are adjacent (there is an edge starting in \code{from} and ending in \code{to}}
\item{isConnected \code{signature(object = "pathway")}}{Returns \code{TRUE} if a pathway contains only one connected component}
\item{isDirected \code{signature(object = "pathway")}}{Returns \code{TRUE} if all edges in a pathway are directed}
\item{edgemode \code{signature(object = "pathway")}}{Returns the type of edges in a pathway: \code{directed}, \code{undirected} or both}
\item{numEdges \code{signature(object = "pathway")}}{Returns the number of edges in a pathway}
\item{numNodes \code{siganture(object = "pathway")}}{Returns the number of nodes in a pathway}
\item{edgeNames \code{signature(object = "pathway")}}{Returns the names of the edges in a following format: starting node ~ ending node}



}
  All of the methods below return an object of class \code{pathway} with modified topology.
\describe{
\item{intersection \code{signature(x = "pathway", y = "pathway")}}{compute the intersection of the two supplied graphs. They must have identical nodes.}

\item{join \code{signature(x = "pathway", y = "pathway")}}{returns the joining of the two graphs. It is similar to \code{intersection} but does not require the identical nodes}

\item{union \code{signature(x = "pathway", y = "pathway")}}{compute the union of the two supplied graphs. They must have identical nodes.}

\item{subGraph \code{signature(snodes = "character", graph = "pathway")}}{Given a set of nodes and a pathway this function creates and returns subgraph with only the supplied nodes and any edges between them  }          

\item{clearNode \code{signature(node = "character", object = "pathway")}}{Clears all edges incoming and outgoing edges from \code{node}}

\item{removeEdge \code{signature(from = "character", to = "character", graph = "character")}}{removes all directed edges starting in \code{from} and ending in \code{to} and undirected edges between \code{from} and \code{to} }

\item{removeNode \code{signature(node = "character", object = "pathway")}}{removes node(s) \code{node} from a pathway \code{object}}

\item{nodes<- \code{signature(x = "pathway", value = "character")}}{sets node labels of pathway \code{object} to \code{value}}
\item{convertIdentifiers \code{signature(x = "pathway", to = "character")}}{converts the node identifiers/labels in a pathway. \code{to} is the name of one of the columns provided by an Annotation package (e.g.\code{"SYMBOL"}}

}


}

\value{
A \code{pathway}
}


\keyword{classes}
