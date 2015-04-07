\name{TopologicalAnalysis}
\docType{methods}
\alias{TopologicalAnalysis}
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



\title{
Methods for basic topological analysis of \code{pathway} class from \pkg{graphite}
}
\description{Describes the methods for basic topological analysis of \code{pathway}. They were defined as generic in \code{graph} for \code{graph} class. They were implemented for \code{pathway} in this package}
\section{Methods}{
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
}

\author{
Ivana Ihnatova
}

\keyword{ methods }

