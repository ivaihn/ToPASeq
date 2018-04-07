.preparePerms <- function(de, all, nperm=1000, method) 
{
    if(method == "PRS") 
    {
        ind <- as.integer(all %in% names(de))
        perms.ind <- replicate(nperm, sample(ind))
        rownames(perms.ind) <- all
        perms <- apply(perms.ind, 2, function(x) {
            x[x == 1] <- sample(de)
            x
        })
    }
    #else if (method == "PWEA") perms <- replicate(nperm, test(x, sample(group))$stat)
    return(perms)
}

.preparePathways <- function(pwys, method, both.directions=TRUE, 
    genes, maxNodes=NULL, minEdges=5, commonTh=5) 
{
    stopifnot(is.list(pwys) || is(pwys, "PathwayList")) 
    
    N <- length(pwys)
    if(!is.null(maxNodes)) pwys <- .bigPaths(pwys, maxNodes)
    if(!is.null(minEdges)) pwys <- .fewEdges(pwys, minEdges)
    
    pwys <- lapply(pwys, function(p) .transformPathway(p, method, both.directions))
    if(!is.null(commonTh)) pwys <- .commonGenes(pwys, genes, commonTh)
    
    nr.filtered <- N - length(pwys)
    if(nr.filtered) message(paste(nr.filtered, "pwys were filtered out"))
    return(pwys)
}

.transformPathway <- function(pwy, method, both.directions=TRUE) 
{
    stopifnot(class(pwy) == "Pathway")
    
    if (method == "TAPPA") 
    {
        pwy <- as(pwy, "graphNEL")
        pwy <- as(pwy, "matrix")
        pwy <- pwy + t(pwy)
        pwy[pwy > 1] <- 1
        pwy[lower.tri(pwy)] <- 0
        diag(pwy) <- 1
    }
    else if(method %in% c("PRS", "PWEA")) 
    {
        pwy <- .buildGraphNEL(nodes(pwy), edges(pwy), both.directions)
        if(method == "PRS") pwy <- as(pwy, "matrix")
    }
    return(pwy)
}

.symmetricEdges <- function(m) 
{
    rind <- (m[, 3] == "undirected") & (m[, 1] != m[, 2])
    undirected <- m[rind, 2:1, drop = FALSE]
    full <- m[, -3, drop = FALSE]
    if(nrow(undirected)) 
    {
        dimnames(full) <- NULL
        full <- rbind(full, undirected)
    } 
    return(full)
}

.edLi <- function(n, e) 
    sapply(n, function(n) list(edges = e[e[, 1] == n, 2]), 
        simplify=FALSE, USE.NAMES = TRUE)


.buildGraphNEL <- function(n, e, both.directions=TRUE) 
{
    if (nrow(e) == 0) g <- new("graphNEL", n, list(), "directed") 
    else 
    {
        e <- as.matrix(e[,c("src", "dest", "direction")])
        n <- vapply(n, function(s) unlist(strsplit(s, ":"))[2], 
                                        character(1), USE.NAMES=FALSE)
        if(both.directions) e <- .symmetricEdges(e)
        g <- new("graphNEL", n, .edLi(n, e), "directed")
        graph::edgeDataDefaults(g, "edgeType") <- "undefined"
        graph::edgeData(g, e[, 1], e[, 2], "edgeType") <- "directed"
    }
    return(g)
}

.fewEdges <- function(pwys, minEdges) 
    Filter(function(p) nrow(edges(p)) > minEdges, pwys)

.bigPaths <- function(pwys, maxNodes)
    Filter(function(p) length(nodes(p)) <= maxNodes, pwys)

.commonGenes <- function(pwys, genes, threshold)
    Filter(function(p) length(intersect(rownames(p), genes)) >= threshold, pwys)

#.dagOnly <- function(pwys) Filter(function(p) gRbase::is.DAG(p), pwys)

