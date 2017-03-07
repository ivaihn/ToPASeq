parseKEGG<-function(pathfolder, method){
 paths<-list.files(pathfolder)
 paths<-paths[regexpr(".xml", paths)!=-1]
 message(paste("Found", length(paths), "KGML files"))
 paths<-paste0(pathfolder,paths)
 
 pt<-list()
 if (method=="SPIA") {
   pt<-lapply(paths, function(x) parseAsSPIA(x))
   names(pt)<-substr(paths, nchar(paths)-8, nchar(paths)-4)
   pt<-Filter(function(p) !is.null(p), pt)} else
 if (method=="DEGraph") {pt<-parseAsDEGraph(paths)} else
 if (method=="PRS") {
   pt<-lapply(paths, function(x) {return(parseAsPRS(x))})
   names(pt)<-substr(paths, nchar(paths)-8, nchar(paths)-4)
} else
 if (method=="CePa") {
   pt<-lapply(paths, function(x) {parseAsCePa(x)})
   pt<-lapply(pt, function(x) as(x,"graphNEL"))
   pt<-lapply(seq_len(length(pt)), function(i) list(pt[[i]], substr(paths[i], nchar(paths[i])-8, nchar(paths[i])-4)))
   names(pt)<-substr(paths, nchar(paths)-8, nchar(paths)-4)
   pt <- Filter(function(p) numEdges(p[[1]])>0, pt)
 } else
 if (method=="clipper") {
 pt<-lapply(paths, function(f) {return(KEGG2Pathway(f, nongene="propagate"))})
 names(pt)<-substr(paths, nchar(paths)-8, nchar(paths)-4)
 pt <- Filter(function(p) !is.null(p), pt)
 pt<-lapply(pt, function(f) {
  nodes(f)<-gsub("hsa:","", nodes(f))
  f@identifier<-"ENTREZ"
  return(f)
 })
 } else
  
 stop("Method: ", method," is not supported")
 return(pt)
}


#SPIA
parseAsSPIA<-function(file) {
LT<-list()
kgml<-parseKGML(file)
nodes<-nodes(kgml)

nodes<-lapply(nodes, function(x) {
 if (x@type=="gene") {
  out<-unlist(lapply(x@name, function(y) strsplit(y,":")[[1]][2]))
} else
 if (x@type=="group") {
   out<-sapply(nodes[x@component], function(y) y@name)
   out<-unlist(out)
   out<-unlist(lapply(out, function(y) strsplit(y,":")[[1]][2]))
 } else out<-c()
return(out)
})

edges<-edges(kgml)
if (length(edges)>0) {
edges<-lapply(edges, function(x) {
st <- sapply(x@subtype, function(y) y@name)
if (length(st)>0) {
st <- sort(st)
st <- st[!st %in% c("", " ")]
st <- paste(st, collapse = "_")
if (st == "dephosphorylation_inhibition") {st = "inhibition_dephosphorylation"}
if (st == "indirect effect_inhibition") {st = "inhibition_indirect effect"}
if (st == "activation_indirect effect_phosphorylation") {
    st = c("activation_indirect effect", "activation_phosphorylation")
    } 
} else st<-""
return(expand.grid(from=x@entry1ID, to=x@entry2ID, type=st))
})
edges<-Reduce(rbind, edges)
edges[,3]<-as.character(edges[,3])

nod<-unique(unname(unlist(unname(nodes))))
LT<-lapply(unique(edges[,3]), function(et) {
  m<-matrix(0, nrow=length(nod), ncol=length(nod), dimnames=list(nod, nod))
  edgessub<-edges[edges[,3]==et,, drop=FALSE]
  for (i in seq_len(nrow(edgessub))) m[nodes[[edgessub[i,1]]],nodes[[edgessub[i,2] ]] ]<-1
  return(m)
})
names(LT)<-unique(edges[,3])
}
LT[["nodes"]] <- unname(unlist(unname(nodes)))
LT[["title"]] <- kgml@pathwayInfo@title
LT[["NumberOfReactions"]] <- length(kgml@reactions)
                    
rel <- c("activation", "compound", "binding/association", 
            "expression", "inhibition", "activation_phosphorylation", 
            "phosphorylation", "inhibition_phosphorylation", 
            "inhibition_dephosphorylation", "dissociation", "dephosphorylation", 
            "activation_dephosphorylation", "state change", "activation_indirect effect", 
            "inhibition_ubiquination", "ubiquination", "expression_indirect effect", 
            "inhibition_indirect effect", "repression", "dissociation_phosphorylation", 
            "indirect effect_phosphorylation", "activation_binding/association", 
            "indirect effect", "activation_compound", "activation_ubiquination")
betas = c(1, 0, 0, 1, -1, 1, 0, -1, -1, 0, 0, 1, 0, 1, -1, 0, 1, -1, -1, 0, 0, 1, 0, 1, 1)
names(betas) = rel
if (LT[["NumberOfReactions"]]>0) return(NULL) else return(LT)
}


SPIA4analysis<-function(pathways, EdgeAttrs=NULL) {
if (is.null(EdgeAttrs)) EdgeAttrs<-makeDefaultEdgeData()
out<-lapply(pathways, function(x)  {
sel<-which(EdgeAttrs$beta$rel %in% names(x))
if (length(sel)==0) res<-NULL else
res<-getdatp(x[EdgeAttrs$beta$rel[sel]], EdgeAttrs$beta$rel[sel], EdgeAttrs$beta$beta[sel])
return(res)
})
message("Pathways: ", paste(names(which(sapply(out, is.null))), collapse=" "), "did not contained any edges of considered type:", paste(EdgeAttrs$beta$rel, collapse=" "),")")
out<-Filter(function(p) !is.null(p), out)
out<-Map(function(x,y) list(x, y), out, names(out))
return(out)
}

SPIA2graph<-function(pathways){
 out<-lapply(pathways, function(x) {
  edgs<-sapply(x, is.matrix) 
  if (sum(edgs)==0) g<-NULL else {
    g<-Reduce("+", x[edgs])
    g<-as(t(g), "graphNEL")
  }
  return(g)
})
out<-Filter(function(p) !is.null(p), out)
}
#DEGraph
parseAsDEGraph<-function(pathnames){
paths_DEGraph <- lapply(pathnames, FUN = function(pathname) {
        pw <- parseKGML(pathname)
        pwInfo <- getPathwayInfo(pw)
        gr <- KEGGpathway2Graph(pw, genesOnly = TRUE, expandGenes = TRUE)
        gr@graphData$info <- pwInfo
        gr@graphData$label <- getTitle(pwInfo)
        gr
    })
pathnames<-substr(pathnames, nchar(pathnames)-8, nchar(pathnames)-4)
paths_DEGraph<-lapply(1:length(paths_DEGraph), function(i) list(paths_DEGraph[[i]], pathnames[i]))
names(paths_DEGraph)<-pathnames
return(paths_DEGraph)
}


#PRS
parseAsPRS<-function(file) {
kgml<-parseKGML(file)
nodA<-sapply(kgml@nodes, KEGGgraph::getName)
nodA<-sapply(nodA, paste, collapse=" ")
nodI<-sapply(kgml@nodes, KEGGgraph::getEntryID)
nodt<-sapply(kgml@nodes, KEGGgraph::getType)

nodTab<-cbind(nodI, nodA, nodt)

nod<-nodA[nodt=="gene"]
nod<-unique(nod)

adjmat<-matrix(0,length(nod), length(nod), dimnames = list(nod,nod))

if (length(kgml@edges)>0  ) {
  edg<-t(sapply(kgml@edges, getEntryID))
  edgN<-apply(edg, 2, function(x) nodTab[match(x,nodTab[,1]),2])
  edgN<-matrix(edgN, ncol=2)
  edgT<-apply(edg, 2, function(x) nodTab[match(x,nodTab[,1]),3])
  edgT<-matrix(edgT, ncol=2)
  edgT<-apply(edgT, 1, paste, collapse="")
  sel<-edgT=="genegene"
  edgF<-edgN[sel,, drop=F]
  edgF<-edgF[!duplicated(edgF),, drop=FALSE]
} else
  edgF<-matrix(ncol=2,nrow=0)
 if (length(kgml@reactions)>0) {
   reactions<-reactionsPRS(kgml)
 } else 
  reactions<-matrix(ncol=2,nrow=0)

if (nrow(reactions)>0) E<-edgF[which(apply(edgF, 1, paste0, collapse="") %in% apply(reactions, 1, paste0, collapse="")),, drop=FALSE] else E<-edgF

#E<-rbind(edgF, reactions)  
 E<-E[regexpr("hsa", E[,1])!=-1 & regexpr("hsa", E[,2])!=-1,, drop=FALSE]
if (nrow(E)>0) 
 for (i in seq_len(nrow(E))) {
  adjmat[E[i,1], E[i,2]]<-1
}
diag(adjmat)<-0

return(list(adjmat, substr(file, 4,8)))
}

reactionsPRS<-function(kgml){
  nodN<-sapply(kgml@nodes, getName)
  nodR<-sapply(kgml@nodes, function(x) x@reaction)
  hasR<-length(kgml@reactions)>0
  
  if (hasR) {
    reaN<-nodN[match(sapply(kgml@reactions, getName), nodR)]
    reaN<-sapply(reaN, paste, collapse=" ")
    reaS<-sapply(kgml@reactions, function(x) x@substrateName)
    reaP<-sapply(kgml@reactions, function(x) x@productName)
    reaT<-sapply(kgml@reactions, function(x) x@type)
    reactions<-Map(expand.grid,reaN, reaS, reaP, reaT)
    reactions<-Reduce( rbind, reactions)
    colnames(reactions)<-c("gene", "substrate", "product", "type")
    reactions<-as.matrix(reactions)
    reactions<-reactions[reactions[,1]!="" ,] #& reactions[,1]!=reactions[,2]
    #   reactions<-cbind(
    #     gene=nodN[match(sapply(kgml@reactions, getName), nodR)],
    #     substrate=sapply(kgml@reactions, function(x) x@substrateName),
    #     product=sapply(kgml@reactions, function(x) x@productName)
    #   )
    reactions<-rbind(reactions,reactions[reactions[,4]=="reversible", c(1, 3, 2, 4)])
    ind<-which(sapply(reactions[,"product"], function(x) reactions[,"substrate"]==x) , arr.ind = TRUE)[,c(2,1), drop=FALSE]
    out<-cbind(reactions[ind[,1],1], reactions[ind[,2],1]) 
    out<-out[!duplicated(out),, drop=FALSE]
  } else out<-matrix(ncol=2, nrow=0)
  return(out)
}

#PWEA
parseAsPWEA<-function(file){
  kgml<-parseKGML(file)
  nodA<-sapply(kgml@nodes, KEGGgraph::getName)
  nodI<-sapply(kgml@nodes, KEGGgraph::getEntryID)
  nodt<-sapply(kgml@nodes, KEGGgraph::getType)
  
  nodTab<-cbind(nodI, nodA, nodt)
  
  adjmat<-matrix(0,length(nodI), length(nodI), dimnames = list(nodI,nodI))
  
  if (length(kgml@edges)>0) {
    edg<-t(sapply(kgml@edges, KEGGgraph::getEntryID))
    edgF<-rbind(edg, edg[,2:1])
    
  } else
    edgF<-matrix(ncol=2,nrow=0)
  if (length(kgml@reactions)>0) {
    reactions<-reactionsPWEA(kgml)
    reactions<-rbind(reactions, reactions[,2:1])
  } else 
    reactions<-matrix(ncol=2,nrow=0)
  
  E<-rbind(edgF, reactions)
  for (i in seq_len(nrow(E))) {
    adjmat[as.character(E[i,1]), as.character(E[i,2])]<-1
  }
  
  return(adjmat)
}

reactionsPWEA<-function(kgml){
  nodN<-sapply(kgml@nodes, KEGGgraph::getName)
  nodR<-sapply(kgml@nodes, function(x) x@reaction)
  hasR<-length(kgml@reactions)>0
  
  if (hasR) {
    
    reaN<-sapply(kgml@reactions, function(x) x@name)
    reaG<-sapply(reaN, function(x) names(nodR[regexpr(x,nodR)!=-1 & !is.na(nodR)]))
    #reaG<-sapply(reaN, function(x) names(nodR[nodR==x & !is.na(nodR)]))
    
    reaS<-lapply(kgml@reactions, function(x) x@substrateName)
    reaP<-lapply(kgml@reactions, function(x) x@productName)
    reaSP<-Map(function(x,y) {
      out<-c(x,y)
      if (length(out)>2) out<-out[1:2]
      return(out)
      }, reaS, reaP)
    
    reactions<-Map(expand.grid,reaG, reaSP)
    reactions<-Reduce( rbind, reactions)
    reactions<-as.matrix(reactions)
    reactions<-reactions[,2:1]
    
    nodN2<-sapply(nodN, paste, collapse=" ")
    reactions[,1]<-names(nodN2)[match(reactions[,1], nodN2)]
    out<-reactions
  } else out<-matrix(ncol=2, nrow=0)
  return(out)
}

#CePa
parseAsCePa<-function(file){
  kgml<-parseKGML(file)
  nodA<-sapply(kgml@nodes, KEGGgraph::getName)
  nodA<-sapply(nodA, paste, collapse=" ")
  
  complex<-sapply(kgml@nodes,KEGGgraph::getComponent)
  nodA<-sapply(complex, function (x) paste(nodA[x], collapse=" ") )
  
  nodA<-nodA[nodA!=""]
  nodTab<-t(t(nodA))
  adjmat<-matrix(0,length(unique(nodA)), length(unique(nodA)), dimnames = list(unique(nodA),unique(nodA)))
  
   if (length(kgml@edges)>0 & length(kgml@reactions)==0) {
     edg<-t(sapply(kgml@edges, getEntryID))
     edg<-apply(edg, 2, function(x) nodTab[match(x,rownames(nodTab)),1])
     edgF<-matrix(edg, ncol=2)
   } else
    edgF<-matrix(ncol=2,nrow=0)
  if (length(kgml@reactions)>0) {
    reactions<-reactionsCePa(kgml)
  } else 
    reactions<-matrix(ncol=2,nrow=0)
  
  E<-rbind(edgF, reactions)
  for (i in seq_len(nrow(E))) {
    if (all(E[i,] %in%  colnames(adjmat))) adjmat[as.character(E[i,1]), as.character(E[i,2])]<-1
  }

  diag(adjmat)<-0
  return(adjmat)
}

reactionsCePa<-function(kgml){
  nodN<-sapply(kgml@nodes, KEGGgraph::getName)
  nodN<-sapply(nodN, paste0, collapse=" ")
  nodR<-sapply(kgml@nodes, function(x) x@reaction)
  hasR<-length(kgml@reactions)>0
  
  if (hasR) {
    
    reaN<-sapply(kgml@reactions, function(x) x@name)
    reaG<-sapply(reaN, function(x) nodN[regexpr(x,nodR)!=-1 & !is.na(nodR)])
    reaS<-lapply(kgml@reactions, function(x) x@substrateName)
    reaP<-lapply(kgml@reactions, function(x) x@productName)
    
    reactions<-Map(function(A,I,O){
      
      if (length(I)==0) out<-expand.grid(A,O) else out<-expand.grid(A,I)
      rbind(out, expand.grid(I,O))
    },reaG, reaS, reaP)
    reactions<-Reduce( rbind, reactions)
    
    reactions<-as.matrix(reactions)
    
    out<-reactions
  } else out<-matrix(ncol=2, nrow=0)
  return(out)
}

