\name{Pathway-method}
\Rdversion{1.1}
\docType{class}
\alias{edges,Pathway,ANY-method}
\alias{nodes,Pathway-method}
\alias{show,Pathway-method}

\alias{intersection,Pathway,Pathway-method}
\alias{join,Pathway,Pathway-method}
\alias{union,Pathway,Pathway-method}
\alias{subGraph,character,Pathway-method}
\alias{degree,Pathway,character-method}
\alias{degree,Pathway,missing-method}
\alias{numNoEdges,Pathway-method}
\alias{mostEdges,Pathway-method}
\alias{acc,Pathway,character-method}
\alias{connComp,Pathway-method}
\alias{edges,Pathway,character-method}
\alias{isAdjacent,Pathway,character,character-method}
\alias{isConnected,Pathway-method}
\alias{isDirected,Pathway-method}
\alias{edgemode,Pathway-method}
\alias{numEdges,Pathway-method}
\alias{numNodes,Pathway-method}
\alias{edgeNames,Pathway-method}

\alias{clearNode,character,Pathway-method}
\alias{removeEdge,character,character,Pathway-method}
\alias{removeNode,character,Pathway-method} 
\alias{addEdge,character,character,Pathway,numeric-method}
\alias{addNode,character,Pathway,list-method}
\alias{nodes<-,Pathway,character-method}
\alias{changeDirection,character,character,Pathway,character-method} 
\alias{changeDirection}
\alias{changeInteraction,character,character,Pathway,character,logical-method}
\alias{changeInteraction}
\alias{convertIdentifiers,Pathway,character-method}

\title{Class \code{"Pathway"}}
\description{
  This class represents a biological pathway.   \code{changeInteraction} and \code{changeDirection} are a new generic function designed for \code{Pathway} class
}

\section{Methods}{
  \describe{
    \item{edges}{\code{signature(object = "Pathway")}:
      retrieves the \code{data.frame} describing the pathway edges.
    }

    \item{nodes}{\code{signature(object = "Pathway")}: retrieves the
      \code{vector} enumerating the identifiers of the pathway nodes.
    }
  }
  The methods below perfom basic topological analysis of a pathway. They were defined as generic in \code{graph} for \code{graph} class. They were implemented for \code{Pathway} in this package
  \describe{

\item{degree \code{signature(object = "Pathway", Nodes = "character")}}{Returns the number of incoming or outgoing edges for nodes in \code{Nodes} }
\item{degree \code{signature(object = "Pathway", Nodes = "missing")}}{Returns the number of incoming or outgoing edges for all nodes in \code{object} }
\item{numNoEdges \code{signature(objGraph = "Pathway")}}{Returns the number of nodes without any edge}
\item{mostEdges \code{signature(objGrapg = "Pathway")}}{Returns the nodes with most edges}
\item{acc \code{signature(object = "Pathway", index = "character")}}{Returns the set of nodes accessible from nodes in \code{index}. The undirected edges are considered as bidirected (directed in both directions)}
\item{connComp \code{signature(object = "Pathway")}}{Returns the connected components present in a pathway. They are returned as list where each slot refers to one component and contains the relevant nodes. The undirected edges are considered as bidirected (directed in both directions)}
\item{edges \code{signature(object = "Pathway", which = "character")}}{Returns the edges relevant to node(s) in \code{which}}
\item{isAdjacent \code{signature(object = "Pathway", from = "character", to = "character")}}{Returns whether nodes in \code{from} and \code{to} are adjacent (there is an edge starting in \code{from} and ending in \code{to}}
\item{isConnected \code{signature(object = "Pathway")}}{Returns \code{TRUE} if a pathway contains only one connected component}
\item{isDirected \code{signature(object = "Pathway")}}{Returns \code{TRUE} if all edges in a pathway are directed}
\item{edgemode \code{signature(object = "Pathway")}}{Returns the type of edges in a pathway: \code{directed}, \code{undirected} or both}
\item{numEdges \code{signature(object = "Pathway")}}{Returns the number of edges in a pathway}
\item{numNodes \code{siganture(object = "Pathway")}}{Returns the number of nodes in a pathway}
\item{edgeNames \code{signature(object = "Pathway")}}{Returns the names of the edges in a following format: starting node ~ ending node}



}
  All of the methods below return an object of class \code{Pathway} with modified topology.
\describe{
\item{intersection \code{signature(x = "Pathway", y = "Pathway")}}{compute the intersection of the two supplied graphs. They must have identical nodes.}

\item{join \code{signature(x = "Pathway", y = "Pathway")}}{returns the joining of the two graphs. It is similar to \code{intersection} but does not require the identical nodes}

\item{union \code{signature(x = "Pathway", y = "Pathway")}}{compute the union of the two supplied graphs. They must have identical nodes.}

\item{subGraph \code{signature(snodes = "character", graph = "Pathway")}}{Given a set of nodes and a pathway this function creates and returns subgraph with only the supplied nodes and any edges between them  }          

\item{clearNode \code{signature(node = "character", object = "Pathway")}}{Clears all edges incoming and outgoing edges from \code{node}}

\item{removeEdge \code{signature(from = "character", to = "character", graph = "Pathway")}}{removes all directed edges starting in \code{from} and ending in \code{to} and undirected edges between \code{from} and \code{to} }

\item{removeNode \code{signature(node = "character", object = "Pathway")}}{removes node(s) \code{node} from a pathway \code{object}}

\item{nodes<- \code{signature(x = "Pathway", value = "character")}}{sets node labels of pathway \code{object} to \code{value}}
\item{convertIdentifiers \code{signature(x = "Pathway", to = "character")}}{converts the node identifiers/labels in a pathway. \code{to} is the name of one of the columns provided by an Annotation package (e.g.\code{"SYMBOL"}}

}


}




\keyword{classes}
