reduceGraph<-function(graph, reduction){

for ( r in 1:length(reduction))  {
newnodes<-c(setdiff(nodes(graph), reduction[[r]]), names(reduction)[r])

e<-edges(graph)[edges(graph)[,1] %in% reduction[[r]] | edges(graph)[,2] %in% reduction[[r]],]


if (!isComplex(graph, reduction[[r]])) {  
 outedg<-tapply(e[,2], e[,1], function(x) x)[reduction[[r]]]
 inedg<-tapply(e[,1], e[,2], function(x) x)[reduction[[r]]]

 outok<-sapply(1:(length(outedg)-1), function(i) all(outedg[[i]]==outedg[[i+1]]))
 inok<-sapply(1:(length(inedg)-1), function(i) all(inedg[[i]]==inedg[[i+1]]))
 if (!any(outok) | !any(inok)) {
  if (!any(outok)) {cat("Ougoing edges:\n"); print(outedg); cat("\n")}
  if (!any(inok)) {cat("Incoming edges:\n"); print(inedg); cat("\n")}
  stop(paste("Nodes ",paste(reduction[[r]], collapse=", ") ,"can not be joined\n"))
  }
} 

newedges<-e
newedges[,1][newedges[,1] %in% reduction[[r]]]<-names(reduction)[r]
newedges[,2][newedges[,2] %in% reduction[[r]]]<-names(reduction)[r]
newedges<-newedges[!duplicated(newedges),] 
newedges<-newedges[newedges[,1]!=newedges[,2],]
newedges<-rbind(edges(graph)[!(edges(graph)[,1] %in% reduction[[r]] | edges(graph)[,2] %in% reduction[[r]]),],newedges)
rownames(newedges)<-1:nrow(newedges)

graph@nodes<-newnodes
graph@edges<-newedges
}
return(graph)
}
