downstream.nodes <-
function(x, set, de, both.directions){
ed<-graphite::edges(set)
if (both.directions) {
ed<-rbind(ed,ed[ed[,3]=="undirected",c(2,1,3,4)]) 
ed[,3]<-"directed"
}
src<-x
downstream<-ed$dest[ed$src %in% src & ed$direction=="directed" & ed$dest %in% de]

repeat{
src<-downstream
downstream.new<-unique(c(downstream, ed$dest[ed$src %in% src & ed$direction=="directed" & ed$dest %in% de]))
if (length(downstream.new) == length(downstream)) break
downstream<-downstream.new
}
return(downstream)
}

#downstream.nodes<-function(x, set, de){
#set.ig<-igraph.from.graphNEL(pathwayGraph(set))
#subg<-subcomponent(set.ig, x, "out")[-1]
#genes<-V(set.ig)[subg]$name
#genes<-genes[genes %in% de]
#return(genes)
#}

