AdjacencyMatrix2Pathway<-function(adjmat, name="pathway", 
 ident="unknown", database="unknown", species="unknown", date=NULL){
title<-name
nodes<-rownames(adjmat)
ident<-ident
database<-database
species<-species
if (is.null(date)) timestamp<-Sys.Date() else timestamp<-date

adjmat2<-adjmat 

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

e<-data.frame(src=edg[,1], dest=edg[,2], direction=dir, type=type) 

e[,3]<-factor(as.character(e[,3]), levels=c("directed", "undirected"))
e[,4]<-factor(as.character(e[,4]), exclude=NULL)

rownames(e)<-NULL
e[e[,4]=="binding" | e[,4]=="process(indirect)" | e[,4]=="process" ,3]<-"undirected"


gr<-new("pathway", title=title, nodes=nodes, edges=e, 
ident=ident, database=database, species=species, 
timestamp=timestamp)
return(gr)
}


