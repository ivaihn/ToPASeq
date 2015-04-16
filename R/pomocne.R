print.topResult<-function(x,...){
cat("$res\n")
print(head(x$res))
cat("... print truncated\n")
cat("$topo.sig\n")
print(x$topo.sig[1:3])
cat("... print truncated\n")
cat("$degtable\n")
print(head(x$degtable))
cat("... print truncated\n")
}



#getBalance.Account <- function(account) {
#  account$balance
#}
#
#getBalance <- function(object, ...)
#  UseMethod("getBalance")

res.topResult<-function(object){
 return(object$res)
 }
res<-function(object){
UseMethod("res")
}
topo.sig.topResult<-function(object){
 return(object$topo.sig)
 }
topo.sig<-function(object){
UseMethod("topo.sig")
}

degtable.topResult<-function(object){
 return(object$degtable)
 }
degtable<-function(object){
UseMethod("degtable")
}



summary.topResult<-function(object,...)
{
cat("Number of analysed pathways: ")
print(nrow(object$res))
cat("\n")

cat("Begining of the table with the results of pathway analysis: \n")
print(head(object$res))

cat("Topological significance of nodes: \n")
if (is.null(object$topo.sig)) cat("The method does not provides the topological significance of nodes\n") else
sapply(object$topo.sig, function(x) summary(x))

cat("Gene-level statistics: \n")
print(summary(object$degtable))
}
#degtable<-function(object){
#if (any("topResult"==class(object))) return(object$degtable) else
# stop("Object must be from class topResult")
# }



#CommonGenes<-function (pathway, genes, threshold=2)
#{
#    commonNames <- intersect(nodes(pathway), genes)
#    if (length(commonNames) < threshold) {
#        warning("not enough genes in common between pathway \"",
#            pathway@title, "\" and expression data (mismatched identifiers?)")
#        return(TRUE)
#    }
#    else return(FALSE)
#}


filterByTag <-function (tag, l) 
{
    isTagged <- sapply(l, function(x) x[[1]] == tag)
    lapply(l[isTagged], function(x) x[[2]])
}



prepareEdges<-function (m, sym)
{
    ns <- canonicalEdgeNames(m)
    simplified <- matrix(unlist(tapply(1:NROW(m), ns, function(is) mergeEdges(m,
        is))), ncol = 4, byrow = TRUE)
    symmetricEdges(simplified)
   }


mergeEdges<-function (m, is)
{
    h <- m[is[1], ]
    if (length(is) == 1) h
    else {
        if ("undirected" %in% m[is, 3] || any(h[1] != m[is, 1]))
            dir <- "undirected"
        else dir <- "directed"
        c(h[1], h[2], dir, paste(unique(m[is, 4]), collapse = ";"))
    }
}
symmetricEdges<-function (m)
{
    undirected <- m[m[, 3] == "undirected" & m[, 1] != m[, 2],
        c(2, 1, 4), drop = FALSE]
    if (NROW(undirected) > 0) {
        full <- m[, -3, drop = FALSE]
        stopifnot(is.null(dimnames(full)))
        rbind(full, undirected)
    }
    else return(m[, -3, drop = FALSE])
}



canonicalEdgeNames<-function (m)
{
    apply(m, 1, function(e) {
        if (e[1] <= e[2]) paste(e[1], e[2], sep = "|")
        else paste(e[2], e[1], sep = "|")
    })
}

edLi<-function (nodes, edges)
{
    sapply(nodes, function(n) list(edges = edges[edges[, 1] ==
        n, 2]), simplify = FALSE, USE.NAMES = TRUE)
}

buildGraphNEL<-function(nodes, edges, sym)
{
    if (NROW(edges) == 0)
        g <- new("graphNEL", nodes, list(), "directed")
    else {
        edges <- prepareEdges(as.matrix(edges), sym)
        g <- new("graphNEL", nodes, edLi(nodes, edges), "directed")
        graph::edgeDataDefaults(g, "edgeType") <- "undefined"
        edgeData(g, edges[, 1], edges[, 2], "edgeType") <- edges[, 3]
    }
    return(g)
}



#Rgraphviz
getRenderPar<-function (g, name, what = c("nodes", "edges", "graph")) 
{
    what <- match.arg(what)
    nms <- switch(what, nodes = nodes(g), edges = edgeNames(g, 
        recipEdges = graphRenderInfo(g, "recipEdges")), graph = "graph")
    ans <- switch(what, nodes = nodeRenderInfo(g, name), edges = edgeRenderInfo(g, 
        name), graph = graphRenderInfo(g, name))
    if (!is.null(ans) && !any(is.na(ans))) {
        if (!is.null(names(ans))) 
            ans <- ans[nms]
    }    else {
        default <- parRenderInfo(g, what)[[name]][1]
        if (is.null(default)) 
            default <- graph.par.get(what)[[name]][1]
        if (is.null(ans)) {
            ans <- rep(default, length(nms))
        }        else {
            if (!is.null(default)) 
                ans[is.na(ans)] <- default
            ans <- ans[nms]
        }
    }
    ans
}

renderSpline<-function (spline, arrowhead = FALSE, arrowtail = FALSE, len = 1, 
    col = "black", lwd = 1, lty = "solid", bbox, ...) 
{
    mylty <- as.numeric(lty)
    if (!is.na(mylty)) 
        lty <- mylty
    lapply(spline, lines, col = col, lwd = lwd, lty = lty, ...)
    warn <- FALSE
    xyhead <- tail(bezierPoints(spline[[length(spline)]]), 2)
    if (is.function(arrowhead[[1]])) {
        xy <- list(x = xyhead[2, 1], y = xyhead[2, 2])
        try(arrowhead[[1]](xy, col = col, lwd = lwd, lty = lty))
    }    else {
        warn <- drawHead(arrowhead, xyhead, bbox, col, lwd, lty, 
            len, out = TRUE)
    }
    xytail <- head(bezierPoints(spline[[length(spline)]]), 2)
    if (is.function(arrowtail[[1]])) {
        xy <- list(x = xytail[1, 1], y = xytail[1, 2])
        try(arrowtail[[1]](xy, col = col, lwd = lwd, lty = lty))
    }    else {
        warn <- warn | drawHead(arrowtail, xytail[2:1, ], bbox, 
            col, lwd, lty, len, out = FALSE)
    }
    warn
}
