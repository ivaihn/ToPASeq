AdjacencyMatrix2Pathway<-function(adjmat, name="pathway", 
 ident="unknown", database="unknown", species="unknown", date=NULL){
title<-name
if (is.null(rownames(adjmat)) & is.null(colnames(adjmat)))
 stop("Adjacency matrix must have rownames and/or colnames")
if (!is.null(rownames(adjmat)) & !is.null(colnames(adjmat))) 
 if (!all(rownames(adjmat)==colnames(adjmat))) stop("Rownames and colnames of the matrix do not match")
if (!all( unique(as.numeric(adjmat)) %in% c(0,1,-1) )) stop("The matrix contains values other than 0, 1, -1")
 

ident<-ident
database<-database
species<-species
if (is.null(date)) timestamp<-Sys.Date() else timestamp<-date

adjmat2<-adjmat 
sel<-rowSums(adjmat)>0 | colSums(adjmat)>0
adjmat2<-adjmat[sel, sel]

ids<-which(!vapply(dimnames(adjmat2), is.null, logical(1)))
nodes<-dimnames(adjmat2)[[ids[1]]]

if (isSymmetric(adjmat))  adjmat2[upper.tri(adjmat2, diag=FALSE)]<-0
 
if (any(-1==adjmat2)) {
edg.act<-which(adjmat2==1, arr.ind=TRUE)
edg.act<-cbind(nodes[edg.act[,1]],nodes[edg.act[,2]])
edg.inh<-which(adjmat2==-1, arr.ind=TRUE)
edg.inh<-cbind(nodes[edg.inh[,1]],nodes[edg.inh[,2]])
edg<-rbind(edg.act, edg.inh)
type<-c(rep("process(activation)", nrow(edg.act)), 
        rep("process(inhibition)", nrow(edg.inh)))

} else {
edg<-which(adjmat2==1, arr.ind=TRUE)
edg<-cbind(nodes[edg[,1]],nodes[edg[,2]])
type<-rep("process", nrow(edg))
}

if (isSymmetric(adjmat)) dir<-rep("undirected", nrow(edg)) else
	dir<-rep("directed", nrow(edg))

e<-data.frame(src=edg[,1], dest=edg[,2], direction=dir, type=type, stringsAsFactors=FALSE) 

e[,3]<-factor(as.character(e[,3]), levels=c("directed", "undirected"))
e[,4]<-factor(as.character(e[,4]), exclude=NULL)

rownames(e)<-NULL
e[e[,4]=="binding" | e[,4]=="process(indirect)" | e[,4]=="process" ,3]<-"undirected"


gr<-new("Pathway", id=title, title=title,  edges=e, 
database=database,species=species, identifier=ident,  
timestamp=timestamp)
return(gr)
}

AdjacencyMatrix2pathway<-function(adjmat, name="pathway", 
 ident="unknown", database="unknown", species="unknown", date=NULL) {
 .Deprecated("AdjacencyMatrix2Pathway")
  
 }