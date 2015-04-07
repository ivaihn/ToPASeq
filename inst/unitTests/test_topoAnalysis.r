test_accCpp<-function(){
library(gRbase)

graph <- dag(c("me","ve"),c("me","al"),c("ve","al"),c("al","an"),c("al","st"),c("an","st"),c("a","b","c"), c("c","d"))
AM<-as(graph,"matrix")
res<-accCpp(AM, rownames(AM), c("b")) 
checkEquals(length(res), 1)
checkEquals(res[[1]][1], c(a=1))

res<-accCpp(AM, rownames(AM), c("a")) 
checkEquals(length(res),1)
checkEquals(length(res[[1]]),0)
}

test_estimateCF<-function(){
library(gRbase)

graph <- dag(c("me","ve"),c("me","al"),c("ve","al"),c("al","an"),c("al","st"),c("an","st"),c("a","b1","b2"), c("b2","d"))
graph<-addEdge("d", "b1", graph)
res<- estimateCF(graphNEL2pathway(graph))

checkEquals(length(res),2)
checkEquals(res[[2]][[1]],c("b1","b2"))

}

test_isComplex<-function(){
uG <- ug(~a:b+b:c+c:d+d:e+e:f+f:a)
uG <- addEdge("c","a",uG)
uG <- addEdge("b","d",uG)
uG <- addEdge("a","d",uG)
gr <- graphNEL2pathway(uG)

gr<-changeInteraction(letters[1:4],c("a","b","c","d"), gr, rep("binding",6), FALSE)
checkTrue(isComplex(c("a","b","c"), gr))
checkEquals(isComplex(gr, c("a","b","me")), FALSE)
}

test_makeDefaultEdgeData<-function(){
att<-makeDefaultEdgeData()
checkEquals(length(att),2)
checkEquals(sapply(att, length),c(graphite2SPIA=52,beta=2))
 checkEquals(sapply(att, nrow,c(graphite2SPIA=26,beta=25))
 checkEquals(names(att), c("graphite2SPIA","beta"))
}

test_KEGG2pathway<-function(){
file<-system.file("extdata", "hsa04610.xml", package = "ToPASeq")
gr<-KEGG2pathway(file, nongene="keep")
checkEquals(numNodes(gr),74)
checkEquals(numEdges(gr),111)
checkEquals(gr@ident, "KEGGnative")

gr<-KEGG2pathway(file, nongene="discard")
checkEquals(numNodes(gr),69)
checkEquals(numEdges(gr),105)
checkEquals(gr@ident, "KEGGnative")
}

test_reduceGraph<-function(){
g<-kegg[["Prolactin signaling pathway"]]
cf<-estimateCF(g)
gr<-reduceGraph(g, cf[[2]])
checkTrue(numNodes(gr) <= numNodes(g))
checkTrue(numEdges(gr) <= numEdges(g))

obs <- tryCatch(reduceGraph(g, list(letters[1:3])), error=conditionMessage)
checkIdentical("Gene identifiers dont match", obs)

}

test_topoAnal<-function(){
graph <- dag(c("me","ve"),c("me","al"),c("ve","al"),c("al","an"),c("al","st"),c("an","st"),c("a","b","c"), c("c","d"))
g<-graphNEL2pathway(graph)
checkEquals(degree(g,"a"), c(a=2))
checkEquals(mostEdges(g), "al")
checkEquals(edgemode(g),"directed")
checkTrue(isDirected(g))
chechTrue(isAdjacent(g,"c","a"))
 checkIdentical(connComp(g), connComp(graph))
}
