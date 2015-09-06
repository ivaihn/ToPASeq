########################
# Topological Analysis #
########################


setMethod("degree",c("Pathway","character"), function(object, Nodes){
if (nrow(object@edges)==0) edgCount=setNames(rep(0, numNodes(object)), nodes(object)) else {
E<-edges(object)
edgCount<-merge(t(t(nodes(object))), table(E[,1]), by.x=1, by.y=1, all=TRUE)
edgCount<-merge(edgCount, table(E[,2]),by=1, all=TRUE)
rownames(edgCount)<-edgCount[,1]
edgCount<-rowSums(edgCount[,-1], na.rm=TRUE)


selfloop<-table(E[,1][E[,1]==E[,2]])
edgCount[names(selfloop)]<-edgCount[names(selfloop)]-1
}

if (!missing(Nodes)>0) edgCount[Nodes] else edgCount
})


setMethod("degree",c("Pathway","missing"), function(object, Nodes){
degree(object, nodes(object))
})



setMethod("numNoEdges", "Pathway", function(objGraph){
sum(degree(objGraph, nodes(objGraph))==0)
})




setMethod("mostEdges", "Pathway", function(objGraph){
names(which.max(degree(objGraph, nodes(objGraph))))
})



setMethod("acc", c("Pathway","character"), function(object, index){
AM<-transformPathway(object, method="PRS")
accCpp(AM, rownames(AM), index)
}
)

setMethod("connComp", "Pathway", function (object){
am<-transformPathway(object, method="PRS")
am<-am+t(am)
am[am>1]<-1
    NL <- nodes(object)
    marked <- rep(0, length(NL))
    names(marked) <- NL
    done <- FALSE
    rval <- vector("list", 1)
    cnode <- 1
    index <- 1
    nused <- numeric(0)
    while (!done) {
        curracc <- names(accCpp(am, NL, NL[cnode])[[1]])
        rval[[index]] <- curracc
        nused <- c(nused, cnode)
        index <- index + 1
        if (length(curracc) > 0) 
            marked[curracc] <- 1
        marked[cnode] <- 1
        cnode <- match(0, marked)
        if (is.na(cnode)) 
            done <- TRUE
    }
    nmR <- NL[nused]
    nc <- length(rval)
    rL <- vector("list", length = nc)
    for (i in seq_len(nc)) rL[[i]] <- c(nmR[[i]], unname(rval[[i]]))
    return(rL)
}
)

setAs("Pathway","graphNEL", function(from) {
out<-buildGraphNEL(nodes(from), edges(from), TRUE)
})

setMethod("edges", c("Pathway","character"), function(object, which){
E<-object@edges
if (!missing(which)) E[E[,1] %in% which | E[,2] %in% which,] else E
})

setMethod("isAdjacent", c("Pathway","character", "character"), function(object,from, to){
E<-edges(object,from)
any(E[,3][E[,2]==to]=="directed", E[,3][E[,1]==to]=="undirected" , E[,3][E[,2]==to]=="undirected" )
})

setMethod("isConnected", "Pathway", function(object){
 length(connComp(object)) == 1
})

#len orientovane hrany, inak FALSE
setMethod("isDirected", "Pathway", function(object){
all(edgemode(object)=="directed")
})

setMethod("edgemode", "Pathway", function(object){
names(which(table(object@edges$direction)>0))
})

setMethod("numEdges", "Pathway",function(object){
nrow(object@edges)})

setMethod("numNodes", "Pathway", function(object){
length(nodes(object))})

setMethod("edgeNames", "Pathway", function(object){
paste(object@edges[,1], object@edges[,2], sep="~")
})

################
# Manipulation #
################

setMethod("intersection", c("Pathway","Pathway"), function(x, y){
if (!all(intersect(nodes(x), nodes(y))==nodes(x))) stop("Supplied graphs must have identical nodes") else {
ex<-x@edges
ey<-y@edges
ex<-apply(ex, 1, function(x) paste(x, collapse=""))
ey<-apply(ey, 1, function(x) paste(x, collapse=""))

w<-which(ex %in% ey)

new("Pathway",   id=paste("Intersect of ", x@id, " and ", y@id, sep=""),
title=paste("Intersect of ", x@title, " and ", y@title, sep=""),
edges=x@edges[w,],
database=mergeDesc(x@database, y@database),      species=mergeDesc(x@species, y@species),
identifier=mergeDesc(x@identifier, y@identifier),   timestamp=Sys.Date())
}})

setMethod("join", c("Pathway", "Pathway"), function(x,y){
u<-rbind(x@edges, y@edges)
u<-u[!duplicated(u),]

new("Pathway",  id=paste("Union of ", x@id, " and ", y@id, sep=""),
title=paste("Union of ", x@title, " and ", y@title, sep=""),
edges=u,
database=mergeDesc(x@database, y@database), species=mergeDesc(x@species, y@species), 
identifier=mergeDesc(x@identifier, y@identifier), timestamp=Sys.Date())
})

setMethod("union",c("Pathway", "Pathway"), function(x,y){
if (!all(intersect(nodes(x), nodes(y))==nodes(x))) stop("Supplied graphs must have identical nodes. Use join() instead.") else {

u<-rbind(x@edges, y@edges)
u<-u[!duplicated(u),]

new("Pathway", id=paste("Union of ", x@id, " and ", y@id, sep=""),
title=paste("Union of ", x@title, " and ", y@title, sep=""),
edges=u,
database=mergeDesc(x@database, y@database), species=mergeDesc(x@species, y@species), 
identifier=mergeDesc(x@identifier, y@identifier), timestamp=Sys.Date())
}
})

setMethod("subGraph", c("character","Pathway"), function(snodes, graph){
if (!all(snodes %in% nodes(graph))) stop("The pathway does not contain the nodes provided")
E<-graph@edges
keepE<-E[,1] %in% snodes & E[,2] %in% snodes
E<-E[keepE,]
rownames(E)<-seq_len(nrow(E))
new("Pathway", id=graph@id, title=graph@title, edges=E,  database=graph@database, species=graph@species,
identifier=graph@identifier,  timestamp=Sys.Date()
)
})

setMethod("clearNode", c("character","Pathway"), function(node, object){
if (!all(node %in% nodes(object))) stop("The pathway does not contain the nodes to be cleared")
E<-object@edges
remE<-E[,1] %in% node | E[,2] %in% node
E<-E[!remE,]
rownames(E)<-seq_len(nrow(E))
new("Pathway", id=object$id, title=object@title, edges=E, 
database=object@database, species=object@species,
identifier=object@identifier,  timestamp=Sys.Date()
)
})
           
setMethod("removeEdge",c("character","character", "Pathway"), function(from, to, graph){
if (!all(c(from,to) %in% nodes(graph))) stop("The pathway does not contain the nodes provided")
E<-graph@edges
remE<-(E[,1] %in% from & E[,2] %in% to & E[,3]=="directed") |
      ( (E[,1] %in% from & E[,2] %in% to) | (E[,1] %in% to & E[,2] %in% from)  & E[,3]=="undirected")

E<-E[!remE,]
rownames(E)<-seq_len(nrow(E))
new("Pathway", id=graph@id, title=graph@title,  edges=E, 
database=graph@database, species=graph@species, 
identifier=graph@identifier, timestamp=Sys.Date()
)
})

setMethod("removeNode", c("character","Pathway"), function(node, object){
if (!all(node %in% nodes(object))) stop("The pathway does not contain the nodes to be removed")
E<-object@edges
remE<-E[,1] %in% node | E[,2] %in% node
E<-E[!remE,]
rownames(E)<-seq_len(nrow(E))
new("Pathway", id=object@id,title=object@title,  edges=E, 
database=object@database, species=object@species,
identifiers=object@identifiers,  timestamp=Sys.Date()
)
})


setMethod("addEdge", c("character","character", "Pathway", "numeric"), function(from, to, graph, weights=1){
extraNodes<-checkNodes(c(from,to), graph)
if (length(extraNodes)>0) stop(paste("The pathway does not contain these nodes", paste(extraNodes, collapse=", ")))

E<-data.frame(src=from, dest=to, direction="unspecified", type="unspecified")
E<-rbind(graph@edges, E)

rownames(E)<-seq_len(nrow(E))
new("Pathway", id=graph@id, title=graph@title,  edges=E, 
 database=graph@database, species=graph@species,
identifier=graph@identifier, timestamp=Sys.Date()
)

})


setMethod("addNode", c("character", "Pathway", "list"), function(node, object, edges){
if (missing(edges)) E<-object@edges else {
if (is.null(names(edges))) stop("edges must be a named list")
nE<-sapply(edges, length)
src<-unlist(lapply(seq_len(length(edges)), function(i) rep(names(edges)[i], nE[i])))
E<-data.frame(src=src, dest=unlist(edges), direction="unspecified", type="unspecified")
E<-rbind(object@edges,E)
rownames(E)<-seq_len(nrow(E))
}
new("Pathway", id=object@id,title=object@title,  edges=E, 
database=object@database, species=object@species,
identifiers=object@identifiers,  timestamp=Sys.Date()
)

})

setGeneric("changeDirection",  function(from, to, graph, type)
    standardGeneric("changeDirection"))
setMethod("changeDirection", c("character", "character", "Pathway", "character"), function(from, to, graph, type){
extraNodes<-checkNodes(c(from,to), graph)
if (length(extraNodes)>0) stop(paste("The pathway does not contain these nodes", paste(extraNodes, collapse=", ")))

E<-graph@edges
E[E[,1] %in% from & E[,2] %in% to,3]<-type
new("Pathway", id=graph@id, title=graph@title, edges=E, 
database=graph@database, species=graph@species,
identifier=graph@identifier,  timestamp=Sys.Date()
)
})

setGeneric("changeInteraction",  function(from, to, graph, interaction, verbose)
    standardGeneric("changeInteraction"))
setMethod("changeInteraction", c("character", "character", "Pathway", "character", "logical"), function(from, to, graph, interaction, verbose=FALSE){
extraNodes<-checkNodes(c(from,to), graph)
if (length(extraNodes)>0) stop(paste("The pathway does not contain these nodes", paste(extraNodes, collapse=", ")))

defaultEdgeAttrs<-makeDefaultEdgeData()
checkInt<-levels(factor(interaction)) %in% defaultEdgeAttrs$graphite2SPIA[,1]
if (!all(checkInt)) warning(paste("Found new interaction type:", paste(levels(factor(interaction))[!checkInt], collapse=", ")))


E<-graph@edges
select<-E[,1] %in% from & E[,2] %in% to
if (verbose) {
cat("Interactions...\n")
E[select,4]
}
E[,4]<-as.character(E[,4])
#if (nrow(E[select,])!= length(interaction)) stop("the number of relevat edges is not equal to the number of interaction type")
E[select,4]<-interaction

if (verbose) {
cat("Were modified as follows...\n")
E[select,4]
}
E[,4]<-factor(E[,4])
new("Pathway",id=graph@id, title=graph@title, edges=E, 
database=graph@database, species=graph@species,
identifier=graph@identifier,  timestamp=Sys.Date()
)

})

setReplaceMethod(f="nodes", signature=c("Pathway","character"),
definition=function(object,value){
if (length(nodes(object))!=length(value)) stop("length of 'value' not equal to number of nodes")
conv.table<-setNames( value, nodes(object))
object<-convertIdentifiersByVector(object, conv.table, "unknown")
return(object)
} )


#######
# AUX #
#######
mergeDesc<-function(x,y){
if (x==y) out<-x else
out<-paste(x,y, sep=", ")
return(out)
}


checkNodes<-function(nodes, object){
extraNodes<-nodes[!nodes %in% nodes(object)]
return(extraNodes)
}


#############


eliminateNode<-function(object, node){
E<-object@edges

inv.edg<-E[,1] %in% node | E[,2] %in% node

E.tmp<-E[inv.edg,]

src<-as.matrix(E.tmp[!E.tmp[,1] %in% node,])
des<-as.matrix(E.tmp[!E.tmp[,2] %in% node,])

combineRows<-function(x,y) {
if (x[3]==y[3] & x[4]==y[4]) out<-c(x[1], y[2], x[3], x[4])
else stop("Incoming and outgoing edge differs in type")
return(out)}

out<-matrix(NA, ncol=4, nrow=nrow(src)*nrow(des))
for (i in seq_len(nrow(src))) {
for (j in seq_len(nrow(des))) {
out[(j-1)*nrow(src)+i,]<-combineRows(src[i,], des[j,])
}}
colnames(out)<-colnames(E)

E<-rbind(E[!inv.edg,], out)

path<-new("Pathway", id=object@id, title=object@title, 
edges=E , database=object@database, species=object@species,
identifier=object@identifier, timestamp=Sys.Date())
return(path)
}

