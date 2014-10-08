
graphNEL2Pathway<-function(graph, name="pathway", 
 ident="unknown", database="unknown", species="unknown", date=NULL){
title<-name
nodes<-nodes(graph)
ident<-ident
database<-database
species<-species
if (is.null(date)) timestamp<-Sys.Date() else timestamp<-date

adjmat<-as(graph,"matrix")
if (!isDirected(graph)) adjmat[upper.tri(adjmat)]<-0

edg<-which(adjmat==1, arr.ind=TRUE)
edg<-cbind(nodes[edg[,1]],nodes[edg[,2]])
nedges<-nrow(edg)


if ("edgeType" %in%  names(graph@edgeData@defaults) )   {
	type<- edgeData(graph, edg[,1], edg[,2],"edgeType") 
  type.split<-sapply(type, function(x) strsplit(x,";", fixed=TRUE))
  type.count<-sapply(type, function(x) length(strsplit(x,";", fixed=TRUE)[[1]]))
  src<-unlist(sapply(1:nedges, function(i) rep(edg[,1][i], each=type.count[i])))
  dest<-unlist(sapply(1:nedges, function(i) rep(edg[,2][i], each=type.count[i])))
  dir<-rep(graph@graphData$edgemode, length(src))
  e<-data.frame(src=src, dest=dest, direction=dir, type= unname(unlist(type.split)), stringsAsFactors=FALSE) 
  m<-match(paste(e[,1], e[,2]),paste(e[,2], e[,1]))
  indpairs<-cbind(1:length(m), m)
  valid<-rep(NA, length(m))
  for (i in 1:nrow(indpairs)) {
  if (!is.na(valid[i])) next else
  if (is.na(valid[i])) valid[i]<-TRUE
  valid[which(indpairs[,1]==indpairs[i,2] & indpairs[,2]==indpairs[i,1])]<-FALSE
  }
  e<-e[valid,]

  } else    {
  dir<-rep(graph@graphData$edgemode, nrow(edg))
	type<-rep("process(indirect effect)", nrow(edg))

e<-data.frame(src=edg[,1], dest=edg[,2], direction=dir, type=type, stringsAsFactors=FALSE) 
}

e[,3]<-factor(as.character(e[,3]), levels=c("directed", "undirected"))
e[,4]<-factor(as.character(e[,4]), exclude=NULL)

rownames(e)<-NULL
e[e[,4]=="binding" | e[,4]=="process(indirect)" | e[,4]=="process" ,3]<-"undirected"



gr<-new("pathway", title=title, nodes=nodes, edges=e, 
ident=ident, database=database, species=species, 
timestamp=timestamp)
return(gr)
}