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

graphite.SPIA<-function (de, all, pathwaySetName, ...)
{
fakeSystemFile<-function (name, package = NULL, ...)
{
    if (package == "SPIA" && name == "extdata")
        getwd()
    else base::system.file(name, package = package, ...)
}

    optArgs <- list(...)
    if (!is.null(optArgs$organism))
        warning("Ignoring the \"organism\" parameter.")
    optArgs$organism <- pathwaySetName
    spiaEnv <- environment(spia)
    fakeEnv <- new.env(parent = spiaEnv)
    assign("system.file", fakeSystemFile, fakeEnv)
    tryCatch({
        environment(spia) <- fakeEnv
        do.call(spia, c(list(de, all), optArgs))[, c(-2, -12)]
    }, finally = {
        environment(spia) <- spiaEnv
    })
}


CheckNames<-function(pathway, expr){

IDmatchsum<-sapply(pathway, function(x) sum(nodes(x) %in% rownames(expr)))
IDmatchmean<-sapply(pathway, function(x) mean(nodes(x) %in% rownames(expr)))
if (sum(IDmatchsum)==0) stop("Gene labels and node labels do not match. Please, correct your gene identifiers\n")
cat(sum(IDmatchsum),"node labels mapped to the expression data\n")
cat("Average coverage", mean(IDmatchmean,na.rm=TRUE)*100,"%\n")
cat(sum(IDmatchsum==0)," (out of ",length(pathway),") pathways without a mapped node\n", sep="")
}

CommonGenes<-function (pathway, genes, threshold=2)
{
    commonNames <- intersect(nodes(pathway), genes)
    if (length(commonNames) < threshold) {
        warning("not enough genes in common between pathway \"",
            pathway@title, "\" and expression data (mismatched identifiers?)")
        return(TRUE)
    }
    else return(FALSE)
}

filterPathwaysByNodeNum<-function (pathways, maxNodes) 
{
    if (!is.null(maxNodes)) 
        pathways <- Filter(function(p) length(nodes(p)) <= maxNodes, 
            pathways)
    return(pathways)
}


lapplyCapturingErrors<-function (l, f) 
{
    log <- lapply(l, function(x) {
        tryCatch(list("ok", f(x)), error = function(e) list("err", 
            e))
    })
    list(results = Filter(Negate(is.null), filterByTag("ok", 
        log)), errors = sapply(filterByTag("err", log), gettext))
}

filterByTag <-function (tag, l) 
{
    isTagged <- sapply(l, function(x) x[[1]] == tag)
    lapply(l[isTagged], function(x) x[[2]])
}

insufficientCommonGenes<-function (pathway, exprGenes) 
{
    commonNames <- intersect(nodes(pathway), exprGenes)
    if (length(commonNames) < 2) {
        warning("not enough genes in common between pathway \"", 
            pathway@title, "\" and expression data (mismatched identifiers?)")
        return(TRUE)
    }
    else return(FALSE)
}



#spiaEdgeType<-function(){
#type=c("binding", "control(In(ACTIVATION))", "control(In(INHIBITION))", "control(Out(ACTIVATION))",
# "control(Out(INHIBITION))", "control(Out(INHIBITION-COMPETITIVE))", "control(Out(ACTIVATION_UNKMECH))", "control(Out(unknown))",
# "control(indirect)", "process", "process(BiochemicalReaction)", "process(activation)", "process(binding/association)", "process(dephosphorylation)",
# "process(dissociation)", "process(expression)", "process(indirect effect)", "process(indirect)",
# "process(inhibition)", "process(missing interaction)", "process(missing)", "process(phosphorylation)",
# "process(repression)", "process(ubiquitination)", "process(methylation)", "process(state change)" )
#spiaType<-c("binding/association", "activation", "inhibition", "activation",
#"inhibition", "inhibition", "activation", "indirect effect",
#"indirect effect", "activation", "activation", "activation",
#"binding/association", "dephosphorylation", "dissociation", "expression",
#"indirect effect", "indirect effect", "inhibition", "indirect effect",
#"indirect effect", "phosphorylation", "inhibition", "ubiquination","inhibition", "ubiquination")
#return(cbind(type=type, spiaType=spiaType))
#}
#graphite:::spiaEdgeType



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