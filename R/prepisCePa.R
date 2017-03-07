## CePa finalna moja implementacia
# centrality, g=graphNEL
indegree<-function(px){
  return(graph::degree(px)$inDegree)
}

outdegree<-function(px){
  return(graph::degree(px)$outDegree)
}

inreach<-function(px){
  sp<-RBGL::floyd.warshall.all.pairs.sp(px)
  sp[is.infinite(sp)]<-NA
  out<-apply(sp, 2, max, na.rm=TRUE) #inreach
  return(out)
}

outreach<-function(px){
  sp<-RBGL::floyd.warshall.all.pairs.sp(px)
  sp[is.infinite(sp)]<-NA
  out<-apply(sp, 1, max, na.rm=TRUE) #outreach
  return(out)
}

between<-function(px){
  return(betweenness(as(px,"matrix"), nodes(px)))
}

CePaWeights<-function(x, method){
  out<-list()
  length(out)<-length(method)
  names(out)<-method
  if (any(method == "equal.weight")) {
    out[["equal.weight"]]<-setNames(rep(1, length(nodes(x))), nodes(x))
  } 
  if (any(method == "in.degree")) {
    out[["in.degree"]]<-indegree(x)
  } 
  if (any(method == "out.degree")) {
    out[["out.degree"]]<-outdegree(x)
  } 
  if (any(method == "degree")) {
    out[["degree"]]<-degree(x, nodes(x))
  } 
  if (any(method == "betweenness")) {
    out[["betweenness"]]<-between(x)
  } 
  if (any(method == "in.reach")) {
    out[["in.reach"]]<-inreach(x)
  } 
  if (any(method == "out.reach")) {
    out[["out.reach"]]<-outreach(x)
  } 
  
  out<-lapply(out, function(x) {
    nz<-x>0
    if (sum(nz)>0) return(x+min(x[nz])/100) else return(x)
  })
  return(out)
}

# priprava permutacii
# globalne
#preparePermsCePa<-function(de, all, nperm){
#  pdiff<-length(de)/length(all)
#  rde<-sapply(1:nperm, function(i) all[as.logical(rbinom(length(all),1, pdiff))])
#}


preparePermsCePap<-function(de, all, nperm, p){
  pdiff<-length(de)/length(all)
  g<-unique(unlist(strsplit(nodes(p), " ", fixed=TRUE)))
  rde<-sapply(1:nperm, function(i) g[as.logical(rbinom(length(g),1, pdiff))])
}


# analyza jednej drahy
CePaSingle<-function(de, all, weights, rde){
  nperm<-length(rde)
  nod<-names(weights[[1]])
  d<-sapply(nod, function(x) {
    g<-strsplit(x, " ", fixed=T)
    #g<-substr(g, 5, nchar(g))
    as.numeric(any(g %in% de))
  })
  #d<-setNames(as.numeric(nod %in% de), nod)
  
  ps<-sapply(weights, function(w) sum(w*d))
  
  rps<-sapply(rde, function(x) {
    dr<-sapply(nod, function(y) {
      g<-strsplit(y, " ", fixed=TRUE)
      #g<-substr(g, 5, nchar(g))
      as.numeric(any(g %in% x))
    })

    #dr<-setNames(as.numeric(nod %in% x), nod)
    psr<-sapply(weights, function(w) sum(w*dr))
    return(psr)
  })
  #psr<<-rps
  p<-rowSums(rps>=matrix(ps, nrow=nrow(rps), ncol=ncol(rps)))/nperm
  return(p)
}

#kegg2<-lapply(kegg, pathwayGraph)
#"betweenness",
CePaORA<-function(de, all, paths, nperm) {
  cen<-c("equal.weight","in.degree","out.degree",
         "in.reach", "out.reach", "betweenness")  
CenList<-lapply(paths, function(x) {CePaWeights(x[[1]], cen)})
#rde<-preparePermsCePa(as.character(de), as.character(all), nperm)
#moje<- t(sapply(CenList, function(x) CePaSingle(de, all, x, rde)))

RDE<-lapply(1:length(paths), function(i) preparePermsCePap(as.character(de), as.character(all), nperm, paths[[i]][[1]]))
OUT<- t(sapply(seq_len(length(RDE)), function(i) CePaSingle(de, all, CenList[[i]], RDE[[i]])))
rownames(OUT)<-sapply(paths, function(x) x[[2]])
#return(list(globalPerms=moje, pathwayPerms=MOJE))
return(OUT)
}