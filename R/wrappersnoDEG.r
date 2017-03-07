# wrappers
TAPPAnoDEG<-function(x, group, pathways, type, preparePaths=TRUE, norm.method=NULL, test.method=NULL, test=t.test, normalize=TRUE, verbose=FALSE, both.directions=TRUE, 
   maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL){
gedm<-prepareData(x, group, type, method="TAPPA", norm.method)
if (preparePaths) paths<-preparePathways(pathways, method="TAPPA", both.directions, rownames(gedm[[1]]), maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy ) else paths<-pathways
res<-tappa(gedm[[1]], gedm[[2]], paths, test, normalize, verbose)

out<-list(res=res, topo.sig=NULL, degtable=list())
class(out)<-c(class(out), "topResultE","topResult")
return(out)
}

clippernoDEG<-function(x, group, pathways, type, preparePaths=TRUE, norm.method=NULL, test.method=NULL,  method="mean", testCliques=FALSE, nperm=1000, alphaV=0.05, b=1000, permute=TRUE, portions=FALSE,
   both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL){
gedm<-prepareData(x, group, type, method="clipper", norm.method)
if (preparePaths) paths<-preparePathways(pathways, method="clipper", both.directions, rownames(gedm[[1]]), maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy ) else paths<-pathways
 res<-CLIPPER(paths, gedm[[1]], gedm[[2]], method, testCliques, nperm, alphaV, b, permute, portions)

out<-list(res=res$results[1:2], topo.sig=res$results[[3]], degtable=list())
class(out)<-c(class(out), "topResultC","topResult")
return(out)
}

TopologyGSAnoDEG<-function(x, group, pathways, type, preparePaths=TRUE, norm.method=NULL, test.method=NULL , method="mean", alpha=0.05, testCliques=FALSE, ..., 
   both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL ){
gedm<-prepareData(x, group, type, method="TopologyGSA", norm.method)
if (preparePaths) paths<-preparePathways(pathways, method="TopologyGSA", both.directions,rownames(gedm[[1]]), maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy ) else paths<-pathways
res<-topologyGSA(gedm[[1]], gedm[[2]], paths, method, alpha, ... )
if (testCliques) {
message("Testing cliques:\n")
 topo.sig<-testCliques(gedm[[1]], gedm[[2]], paths, method, alpha, ...)
 } else topo.sig<-NULL

out<-list(res=res, topo.sig=topo.sig, degtable=list())
class(out)<-c(class(out), "topResultC","topResult")
return(out)
}

DEGraphnoDEG<-function(x, group, pathways, type, preparePaths=TRUE, norm.method=NULL, test.method=NULL, overall="biggest", useInteractionSigns=TRUE, EdgeAttrs=NULL,
   both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL){
gedm<-prepareData(x, group, type, method="DEGraph", norm.method)

if (preparePaths) {
if (useInteractionSigns) { 
   if (is.null(EdgeAttrs)) EdgeAttrs<-makeDefaultEdgeData()
   paths<-preparePathways(pathways, method="DEGraph", both.directions,rownames(gedm[[1]]), maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy, EdgeAttrs)
   } else
   paths<-preparePathways(pathways, method="DEGraphNoSigns", both.directions, rownames(gedm[[1]]),maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy)
} else paths<-pathways
   
res<-degraph(gedm[[1]], gedm[[2]], paths, overall)

out<-list(res=res, topo.sig=NULL, degtable=list())
class(out)<-c(class(out), "topResultE","topResult")
return(out)
}

PWEA<-function(x, group, pathways, type, preparePaths=TRUE, norm.method=NULL, test.method=NULL, tif=NULL, alpha=0.05, nperm=1000,  ncores=NULL, 
  both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL){
gedm<-prepareData(x, group, type, method="PWEA", norm.method, test.method, nperm, ncores)
if (preparePaths) paths<-preparePathways(pathways, method="PWEA", both.directions, rownames(gedm[[1]]),maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy ) else paths<-pathways
if (is.null(tif) & type=="DEtable") stop("Argument 'tif' is missing. Calculate TIF with prepareTIF() first.")
if (is.null(tif)) tif<-prepareTIF(paths, x, alpha)
res<-pwea(gedm[[1]], gedm[[2]], tif, paths, alpha)
deg.table<-gedm[[1]]
topo.sig<-lapply(paths, function(p) {
g<-nodes(p[[1]])
g<-g[g %in% gedm[[1]]$ID]
tif[[p[[2]]]][g]

# exprs.valid<-extractsubset(x, p[[1]])$x
# return(TIF(x, exprs.valid))
 }
 )
out<-list(res=res, topo.sig, degtable=deg.table)
class(out)<-c(class(out), "topResultW", "topResult")
return(out)
}

SPIAnoDEG<-function(x, group, pathways, type, preparePaths=TRUE, norm.method=NULL, test.method=NULL, p.th=0.05, logFC.th=2, nperm=1000, combine="fisher", 
 both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL){
degs<-prepareData(x, group, type, method="SPIA", norm.method, test.method, p.th=p.th, logFC.th=logFC.th)

if (preparePaths) paths<-preparePathways(pathways, method="SPIA", both.directions, degs[[2]], maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy ) else paths<-pathways
res<-spia(degs[[1]], degs[[2]], paths, nperm, combine)

topo.sig<-collectWeightsSPIA(degs[[1]], degs[[2]], paths)
out<-list(res=res, topo.sig=topo.sig, degtable=list())
class(out)<-c(class(out), "topResultW","topResult")
return(out)
}

PRSnoDEG<-function(x, group, pathways, type, preparePaths=TRUE, norm.method=NULL, test.method=NULL, p.th=0.05, logFC.th=2, nperm=1000,
 both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL){
degs<-prepareData(x, group, type, method="PRS", norm.method, test.method, p.th=p.th, logFC.th=logFC.th)
if (preparePaths) paths<-preparePathways(pathways, method="PRS", both.directions, degs[[2]], maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy ) else paths<-pathways
res<-prs(abs(degs[[1]]), degs[[2]], paths, nperm)


topo.sig<-collectWeightsPRS(degs[[1]], degs[[2]], paths)
out<-list(res=res, topo.sig=topo.sig, degtable=list())
class(out)<-c(class(out), "topResultW","topResult")
return(out)
}


CePanoDEG<-function(x, group, pathways, type, preparePaths=TRUE, norm.method=NULL, test.method=NULL, p.th=0.05, logFC.th=2, nperm=1000, 
 both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL){
degs<-prepareData(x, group, type, method="CePa", norm.method, test.method, p.th=p.th, logFC.th=logFC.th)
degs[[1]]<-names(degs[[1]])
if (preparePaths) paths<-preparePathways(pathways, method="CePa", both.directions, degs[[2]], maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy ) else paths<-pathways
res<-CePaORA(degs[[1]], degs[[2]], paths, nperm)


 cen<-c("equal.weight","in.degree","out.degree",
         "in.reach", "out.reach", "betweenness")  
topo.sig<-lapply(paths, function(x) {CePaWeights(x[[1]], cen)})

out<-list(res=res, topo.sig=topo.sig, degtable=list())
class(out)<-c(class(out), "topResultW","topResult")
return(out)
}