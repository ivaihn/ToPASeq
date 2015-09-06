estimateCF<-function(graph){
if (!any(class(graph)=="Pathway")) stop("graph must be of class 'Pathway',found:", paste(class(graph)))
#complexex
#graph<-as(graph,"pathway")
e<-edges(graph)
n<-length(nodes(graph))
nod<-nodes(graph)


e.sub<-e[e[,3]=="undirected" & e[,4]=="binding",]
compNod<-unique(c(e.sub[,1],e.sub[,2]))

comp<-getCliques(subGraph(compNod, buildGraphNEL(nodes(graph), edges(graph), TRUE)))
comp[sapply(comp, length)==2]<-NULL
if (length(comp)>0) names(comp)<-paste("complex", seq_len(length(comp)),sep="")


#families

outedg <- tapply(e[, 2], e[, 1], function(x) x)[nodes(graph)]; names(outedg)<-nodes(graph)
inedg <- tapply(e[, 1], e[, 2], function(x) x)[nodes(graph)]; names(inedg)<-nodes(graph)

noout<-sapply(outedg, is.null)
noin<-sapply(inedg, is.null)

outMatch<-outer(outedg, outedg, function(x,y) vapply(seq_along(x), function(i) identical(x[[i]],y[[i]]), logical(1)))
outGroups<-apply(outMatch, 1, function(x) paste(which(x), collapse=","))

inMatch<-outer(inedg, inedg, function(x,y) vapply(seq_along(x), function(i) identical(x[[i]],y[[i]]), logical(1)))
inGroups<-apply(inMatch, 1, function(x) paste(which(x), collapse=","))

io<-cbind(inGroups, outGroups)

splint<-function(x){
paste(intersect(strsplit(x,",")[[1]], strsplit(x,",")[[2]]), collapse=",")
}



famMid<-unique(apply(io, 1, splint))
famMid<-lapply(famMid, function(x) nod[as.numeric(strsplit(x,",")[[1]])])


fam<-famMid

single.fam<-sapply(fam, function(x) length(x)==1)
if (sum(single.fam)==length(fam)) fam<-list() else fam[sapply(fam, function(x) length(x)==1)]<-NULL
if (length(fam)>0) {

#names(fam)<-paste("family",seq_len(length(fam)),sep="")
namesF<-sapply(fam, function(X) paste(suppressWarnings(Reduce(f<-function(x,y){x[x %in% y]},strsplit(X,""))), collapse=""))
w<-which(nchar(namesF)==0)
if (length(w)>0) warning("Some representative names could not be estimated")
if (length(namesF) != length(unique(namesF))) warning("Found duplicities in the representative names")
if (length(w)>0) namesF[w]<-paste("family", seq_len(length(w)), sep="")
names(fam)<-namesF
}

out<-list(complexes=comp, families=fam)
return(out)
}
