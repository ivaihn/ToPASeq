# RUnit testing

test_downstream.nodes<-function(){

set.seed(1237)
V<-letters[1:10]
from <- sample(V,7)
 to <- sample(V,7)
 ft <- cbind(from, to)
 ft <- rbind(ft,c("a","b"))
 g <- ftM2graphNEL(ft)
g<-graphNEL2Pathway(g)


checkEquals(length(ToPASeq:::downstream.nodes("a",g,c("a","b","g"),FALSE)),2)
checkEquals(length(ToPASeq:::downstream.nodes("i",g,c("a","b","g"),FALSE)),0)
checkEquals(length(ToPASeq:::downstream.nodes("w",g,c("a","b","g"),FALSE)),0)
}

################
test_convertIdentifiersByVector<-function(){
g<-kegg[["Asthma"]]
conv<-setNames(paste("gene", 1:length(nodes(g)), sep=""), nodes(g))
gc<-convertIdentifiersByVector(g, conv, "dummy")
dum<-matrix(rnorm(20), 2,5, dimnames=list(c("gene1","gene2")))
checkTrue(is.numeric(runDEGraph(gc, dum, c(1,1,1,2,2))[[1]][[1]]))

g<-kegg[["Asthma"]]
conv<-setNames(paste("gene", 1:length(nodes(g)), sep=""), nodes(g))
checkException(convertIdentifiersByVector(g, conv[1], "dummy"), paste("These pathway nodes are missing in the 'conv.table':",
ns[! ns %in% names(conv.table)]) ) 
 }
###############
test_graphNEL2Pathway<-function(){
g1<-kegg[[1]]
g2<-pathwayGraph(g1)
g3<-graphNEL2Pathway(g2)

checkTrue(length(nodes(g1))==length(nodes(g3)))
checkTrue(nrow(edges(g1))==nrow(edges(g3)))

set.seed(1237)
V<-letters[1:10]
from <- sample(V,7)
 to <- sample(V,7)
 ft <- cbind(from, to)
 ft <- rbind(ft,c("a","b"))
 g <- ftM2graphNEL(ft)
gg<-graphNEL2Pathway(g)

checkTrue(length(nodes(g))==length(nodes(gg)))
checkTrue(length(unlist(edges(g)))==nrow(edges(gg)))
checkTrue(nlevels(edges(gg)[,4])==1)
} 