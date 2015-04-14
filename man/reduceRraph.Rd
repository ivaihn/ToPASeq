\name{reduceGraph}
\alias{reduceGraph}

\title{
Function to reduce the pathway graph 
}
\description{
Function simplifies a pathway graph topology. It merges a user specified nodes into a one. The specified set of nodes must be either a gene family or a protein complex. By a gene family we mean a set of genes with same outgoing or incoming edges. On the other hand, a protein complex is a set of nodes with only undirected binding edges between them and the number of edges is equal to the complex size. 
}
\usage{
reduceGraph(graph, reduction)
}

\arguments{
  \item{graph}{
An object of class pathway, a pathway to be reduced
}
  \item{reduction}{
A named list of reductions to be maded. 
}
}

\value{
A \code{pathway}
}
\author{
Ivana Ihnatova
}
\examples{
pathways<-pathways("hsapiens","kegg")["Prolactin signaling pathway"]
pathways<-convertIdentifiers(pathways[[1]], "symbol")


gr<-as(pathways,"pathway")
red<-list(RAS=c("NRAS","KRAS","HRAS"), SHC=c("SHC1", "SHC4","SHC2","SHC3"))
reduced<-reduceGraph(gr, red)
reduced
par(mfrow=c(1,2))

nA<-list(fillcolor=c(NRAS="red", KRAS="red", HRAS="red", SHC1="green", SHC4="green", SHC2="green", SHC3="green"))
plot(as(gr,"graphNEL"), nodeAttrs=nA, attrs=list(node=list(fontsize=30, height=40)), main="Before")

plot(as(reduced,"graphNEL"), 
 nodeAttrs=list(fillcolor=c(RAS="red", SHC="green")), attrs=list(node=list(fontsize=30, height=40)), main="After")

#this throws an error, "RELA", "FOS","NFKB1" is not correct set of genes
\dontrun{
pathways<-pathways("hsapiens","kegg")["Prolactin signaling pathway"]
pathways<-convertIdentifiers(pathways[[1]], "symbol")


gr<-as(pathways,"pathway")
gr<-convertIdentifiers(kegg[["Prolactin signaling pathway"]],"symbol")
red<-list(RAS=c("NRAS","KRAS","HRAS"), SHC=c("RELA", "FOS","NFKB1"))
reduced<-reduceGraph(gr, red)
}

pathways<-pathways("hsapiens","kegg")["p53 signaling pathway"]
pathways<-convertIdentifiers(pathways[[1]], "symbol")


gr<-as(pathways,"pathway")

reduced<-reduceGraph(gr, list(com=c("CCNE1", "CDK2", "CCNE2")))

par(mfrow=c(1,2))

nA<-list(fillcolor=c(CCNE1="red", CDK2="red", CCNE2="red"))
plot(as(gr,"graphNEL"), nodeAttrs=nA, attrs=list(node=list(fontsize=30, height=40)), main="Before")

plot(as(reduced,"graphNEL"), 
 nodeAttrs=list(fillcolor=c(com="red")), attrs=list(node=list(fontsize=30, height=40)), main="After")

}
\keyword{ manip }

