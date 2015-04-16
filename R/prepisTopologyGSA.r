
#########

PMT<-function (x, group, dag, alpha, perms, variance = TRUE, paired = FALSE) 
{
    group<-as.numeric(factor(group))
    y1<-t(x[,group==1])
    y2<-t(x[,group==2])
    
    vs <- c(substitute(y1), substitute(y2), substitute(dag))
    l <- procParams(y1, y2, dag)
    y1 <- l$y1
    y2 <- l$y2
    y <- rbind(y1, y2)
    y.num <- nrow(y1) + nrow(y2)
    l <- .runPathwayVarTest(y1, y2, l$graph, alpha, variance)
    cli.moral <- l$cli.moral
    var.equal <- l$var.equal
    if (paired) {
        t.obs <- .hotePaired(y1, y2, cli.moral)
        stat.perm <- vector("numeric", perms)
        #* for (i in seq_len(perm.num)) stat.perm[i] <- .hotePaired(y1, y2, cli.moral, perm = TRUE)
        p.value <- sum(stat.perm >= t.obs)/perms
        l <- list(p.value = p.value, cli.moral = cli.moral, graph = l$graph, 
            t.value = t.obs)
    }    else {
        s <- .hote(y1, y2, FALSE, cli.moral)
        y1.num <- nrow(y1)
        stat.perm <- vector("numeric", perms)
        for (i in ncol(perms)) {
            
            y1.perm <- y[perms[,i]==1, ]
            y2.perm <- y[perms[,i]==2, ]
            stat.perm[i] <- .hote(y1.perm, y2.perm, FALSE, cli.moral)$t.obs
        }
        p.value <- sum(stat.perm >= s$t.obs)/perms
        l <- list(t.value = s$t.obs, df.mean = s$df, p.value = p.value, 
            lambda.value = l$lambda.value, df.var = l$df, p.value.var = l$p.value, 
            qchisq.value = l$qchisq.value, var.equal = var.equal, 
            cli.moral = cli.moral, graph = l$graph)
    }
    l$var.names <- vs
    attr(l, "class") <- "pathway.mean.test"
    return(l)
}

procParams<-function (y1, y2, dag) 
{
    if (!is.matrix(y1) && !is.data.frame(y1)) 
        stop("y1 is not a matrix or a data frame.")
    else if (!is.matrix(y2) && !is.data.frame(y2)) 
        stop("y2 is not a matrix or a data frame.")
    else if (ncol(y1) != ncol(y2)) 
        stop("y1 and y2 differ in the number of columns (genes)")
    else if (any(colnames(y1) != colnames(y2))) 
        stop("y1 and y2 differ in the column names (gene names)")
    else if (nrow(y1) < 3) 
        stop("y1 should have at least 3 rows (samples)")
    else if (nrow(y2) < 3) 
        stop("y2 should have at least 3 rows (samples)")
    common <- intersect(colnames(y1), nodes(dag))
    if (length(common) < 3) 
        stop("need at least 3 genes in common between expression and dag")
    y1 <- y1[, common, drop = FALSE]
    y2 <- y2[, common, drop = FALSE]
    dag <- subGraph(common, dag)
    graph <- cli(dag)
    list(y1 = y1, y2 = y2, graph = graph)
}

cli<-function (dag) {
    moral <- moralize(dag)
    tg <- triangulate(moral)
    adj.moral <- as(moral, "matrix")
    cli<-getCliques(moral)
    cli<-lapply(cli, function(x) sort(match(x, nodes(dag))))
    list(adj.moral = adj.moral, cli.moral = cli, cli.tg = getCliques(tg), moral = moral, 
        tg = tg)
}


.hote<-function (y1, y2, exact, cliques = NULL) {
    y1.num <- nrow(y1)
    y2.num <- nrow(y2)
    gene.num <- ncol(y1)
    y1.bar <- colMeans(y1)
    y2.bar <- colMeans(y2)
    y1.s <- cov(y1)
    y2.s <- cov(y2)
    y.diff <- y1.bar - y2.bar
    if (!is.null(cliques)) {
        y1.s <- qpIPF(y1.s, cliques)
        y2.s <- qpIPF(y2.s, cliques)
    }
    k <- y1.num + y2.num - gene.num - 1
    if (exact) {
        s <- ((y1.num - 1) * y1.s + (y2.num - 1) * y2.s)/(y1.num + 
            y2.num - 2)
        t2 <- ((y1.num * y2.num)/(y1.num + y2.num)) * (y.diff %*% 
            solve(s) %*% y.diff)
        t.obs <- as.vector(t2 * k/(gene.num * (y1.num + y2.num - 
            2)))
        alpha.obs <- 1 - pf(t.obs, gene.num, k)
        list(alpha.obs = alpha.obs, t.obs = t.obs, df = c(gene.num, 
            k))
    }     else {
        s <- y1.s/y1.num + y2.s/y2.num
        t2 <- as.vector(((y1.num * y2.num)/(y1.num + y2.num)) * 
            (y.diff %*% solve(s) %*% y.diff))
        list(t.obs = t2, df = c(gene.num, k))
    }
}

.hotePaired<-function (y1, y2, cli.moral, perm = FALSE) {


    y1.num <- nrow(y1)
    y.diff <- y1 - y2
    if (perm) {
        signs <- matrix(sample(c(1, -1), y1.num * ncol(y1), replace = TRUE), 
            nrow = y1.num)
        y.diff <- y.diff * signs
    }
    y.bar <- colMeans(y.diff)
    y.centr <- y.diff - y.bar
    y.s <- qpIPF(cov(y.diff), cli.moral)
    t2 <- as.numeric(y1.num * (t(y.bar) %*% solve(y.s) %*% y.bar))
    p <- ncol(y1)
    np <- y1.num - p
    t2 * np/(p * (y1.num - 1))
}

.runPathwayVarTest<-function (y1, y2, graph, alpha, variance, s1 = NULL, s2 = NULL) 
{
    cliques <- graph$cli.moral
    maxCliqueSize <- max(sapply(cliques, length))
    if (nrow(y1) <= maxCliqueSize) 
        stop("y1 should have more than ", maxCliqueSize, " rows (samples)")     else 
    if (nrow(y2) <= maxCliqueSize) 
        stop("y2 should have more than ", maxCliqueSize, " rows (samples)")
    if (is.null(s1) != is.null(s2)) {
        stop("You must provide both s1 and s2 or neither.")
    } else if (is.null(s1) && is.null(s2)) {
        cov <- .estimateCov(y1, y2)
        s1.hat <- qpIPF(cov$s1, cliques)
        s2.hat <- qpIPF(cov$s2, cliques)
        s.hat <- qpIPF(cov$s, cliques)
    }    else {
        s1.hat <- s1
        s2.hat <- s2
        s.hat <- .estimateCovPool(nrow(y1), nrow(y2), s1, s2)
    }
    s1.det <- det(s1.hat)
    s2.det <- det(s2.hat)
    s.det <- det(s.hat)
    lambda.value <- nrow(y1) * log(s.det/s1.det) + nrow(y2) * 
        log(s.det/s2.det)
    df <- (sum(graph$adj.moral)/2) + ncol(y1)
    qchisq.value <- qchisq(1 - alpha, df)
    p.value <- 1 - pchisq(lambda.value, df)
    var.equal <- p.value <= alpha
    res <- list(lambda.value = lambda.value, p.value = p.value, 
        df = df, qchisq.value = qchisq.value, cli.moral = cliques, 
        var.equal = var.equal)
    if (variance) {
        res$s1 <- s1.hat
        res$s2 <- s2.hat
    }
    res$graph <- graph$moral
    return(res)
}

.estimateCov<-function (y1, y2) {
    y1.num <- nrow(y1)
    y2.num <- nrow(y2)
    gene.num <- ncol(y1)
    y1.mean <- apply(y1, 2, mean)
    names(y1.mean) <- NULL
    y1.scal <- y1 - matrix(rep(y1.mean, y1.num * gene.num), y1.num, 
        gene.num, byrow = T)
    s1 <- cov(y1.scal)
    y2.mean <- apply(y2, 2, mean)
    names(y2.mean) <- NULL
    y2.scal <- y2 - matrix(rep(y2.mean, y2.num * gene.num), y2.num, 
        gene.num, byrow = T)
    s2 <- cov(y2.scal)
    cov <- list(s1 = s1, s2 = s2)
    cov$s <- .estimateCovPool(y1.num, y2.num, s1, s2)
    return(cov)
}

.estimateCovPool<-function (y1.num, y2.num, s1, s2) {
    (s1 * (y1.num - 1) + s2 * (y2.num - 1))/(y1.num + y2.num - 
        2)
}


PVT<-function (x, group, dag, alpha, variance = FALSE, s1 = NULL, s2 = NULL) {
    group<-as.numeric(factor(group))
    y1<-t(x[,group==1])
    y2<-t(x[,group==2])
    
    vs <- c(substitute(y1), substitute(y2), substitute(dag))
    l <- procParams(y1, y2, dag)
    y1 <- l$y1
    y2 <- l$y2
    res <- .runPathwayVarTest(y1, y2, l$graph, alpha, variance, 
        s1, s2)
    ns <- nodes(res$graph)
    res$cli.moral <- lapply(res$cli.moral, function(ixs) ns[ixs])
    res$var.names <- vs
    attr(res, "class") <- "pathway.var.test"
    return(res)
}

CMT<-function (x, group, dag, alpha, perms = 1000, paired = FALSE) 
{
    group<-as.numeric(factor(group))
    y1<-t(x[,group==1])
    y2<-t(x[,group==2])
    
    vs <- c(substitute(y1), substitute(y2), substitute(dag))
    l <- procParams(y1, y2, dag)
    y1 <- l$y1
    y2 <- l$y2
    cli.test <- .runCliqueVarTest(y1, y2, l$graph, alpha)
    check <- cli.test$var.equal
    cliques <- cli.test$cliques
    clique.num <- length(cliques)
    alpha.obs <- vector("numeric", clique.num)
    t.obs <- vector("numeric", clique.num)
    for (i in seq_len(clique.num)) {
        cli <- unlist(cliques[i])
        y1.cli <- y1[, cli, drop = FALSE]
        y2.cli <- y2[, cli, drop = FALSE]
        if (length(cli) != 1) {
            if (paired) {
                y.diff <- y1 - y2
                y.bar <- colMeans(y.diff)
                y.centr <- y.diff - y.bar
                y.num <- nrow(y1)
                y.s <- (t(y.centr) %*% y.centr)/y.num
                t2 <- y.num * (t(y.bar) %*% solve(y.s) %*% y.bar)
                p <- ncol(y1)
                np <- y.num - p
                t.value <- t2 * np/(p * (y.num - 1))
                r <- list(alpha.obs = 1 - pf(t.value, p, np), 
                  t.obs = t.value)
            } else if (check[i]) {
                r <- .mult.test(y1.cli, y2.cli, perms)
            } else {
                r <- .hote(y1.cli, y2.cli, TRUE)
            }
            alpha.obs[i] <- r$alpha.obs
            t.obs[i] <- r$t.obs
        } else {
            r <- t.test(y1.cli, y2.cli, paired = paired)
            alpha.obs[i] <- r$p.value
            t.obs[i] <- r$statistic
        }
    }
    l <- list(t.value = t.obs, p.value = as.numeric(alpha.obs))
    if (!paired) {
        l$lambda.value <- cli.test$lambda.value
        l$p.value.var <- cli.test$p.value
        l$var.equal <- check
    }
    l$cliques <- cliques
    l$graph <- cli.test$graph
    l$var.names <- vs
    attr(l, "class") <- "clique.mean.test"
    return(l)
}

.mult.test<-function (y1, y2, perm.num) 
{
    y1.num <- nrow(y1)
    y <- rbind(y1, y2)
    y.num <- nrow(y)
    t.obs <- .hote(y1, y2, FALSE)$t.obs
    stat.perm <- vector("numeric", perm.num)
    for (i in 1:perm.num) {
        ind <- sample(y.num)
        y1.perm <- y[ind[1:y1.num], ]
        y2.perm <- y[ind[(y1.num + 1):y.num], ]
        stat.perm[i] <- .hote(y1.perm, y2.perm, FALSE)$t.obs
    }
    alpha.obs <- sum(stat.perm >= t.obs)/perm.num
    list(alpha.obs = alpha.obs, t.obs = t.obs)
}

CVT<-function (x, group, dag, alpha) {
    group<-as.numeric(factor(group))
    y1<-t(x[,group==1])
    y2<-t(x[,group==2])
    
    vs <- c(substitute(y1), substitute(y2), substitute(dag))
    p <- procParams(y1, y2, dag)
    l <- .runCliqueVarTest(p$y1, p$y2, p$graph, alpha)
    l$var.names <- vs
    attr(l, "class") <- "clique.var.test"
    return(l)
}

.runCliqueVarTest<-function (y1, y2, graph, alpha) {
    cliques <- graph$cli.tg
    maxCliqueSize <- max(sapply(cliques, length))
    if (nrow(y1) <= maxCliqueSize) 
        stop("y1 should have more than ", maxCliqueSize, " rows (samples)") else
    if (nrow(y2) <= maxCliqueSize) 
        stop("y2 should have more than ", maxCliqueSize, " rows (samples)")
    cov <- .estimateCov(y1, y2)
    clique.num <- length(cliques)
    alpha.obs <- rep(0, clique.num)
    lambda.obs <- rep(0, clique.num)
    check <- rep(FALSE, clique.num)
    for (i in seq_along(cliques)) {
        cli <- unlist(cliques[i])
        p <- length(cli)
        s1.hat <- cov$s1[cli, cli, drop = FALSE]
        s2.hat <- cov$s2[cli, cli, drop = FALSE]
        s.hat <- cov$s[cli, cli, drop = FALSE]
        s1.det <- det(s1.hat)
        s2.det <- det(s2.hat)
        s.det <- det(s.hat)
        lambda.obs[i] <- nrow(y1) * log(s.det/s1.det) + nrow(y2) * 
            log(s.det/s2.det)
        alpha.obs[i] <- 1 - pchisq(lambda.obs[i], p * (p + 1)/2)
        if (alpha.obs[i] <= alpha) 
            check[i] <- TRUE
    }
    list(lambda.value = lambda.obs, p.value = alpha.obs, cliques = cliques, 
        var.equal = check, graph = graph$tg)
}



topologyGSA<-function(x, group, pathways, test="mean", alpha, both.directions, ... ){
 test<-switch(test, var = PVT, mean = PMT, stop("invalid test type: ", test))
 #perms<-preparePerms(group=group, nperm=nperm, method="TopologyGSA")
 message("Analysing pathway:\n")
 out<-catchErr(pathways, function(p) {
 cat(p[[2]],"\n")
  test(x, group, p[[1]], alpha, ...)
})

if (length(out[[1]])>0)  {res<-sapply(out[[1]], function(x) unlist(x[sapply(x, function(y) mode(y) %in% c("numeric", "logical"))])) 
#graphs<-sapply(out[[1]], function(x) (x[sapply(x, function(y) mode(y) %in% c("S4", "list"))]))
out[[1]]<-data.frame(t(res))
if (length(out[[1]]>0)) out[[1]]$q.value<-p.adjust(out[[1]][,"p.value"],"fdr")
}
return(out)
}

testCliques<-function(x, group, pathways, test, alpha,  ...){
xx<-list()
if (test=="mean") {
xx<-lapply(pathways, function(p) {
cat(p[[2]],"\n")
 g<-nodes(p[[1]]) 
 g<-g[g %in% rownames(x)]
 res<-CMT(x[g,], group, p[[1]], alpha, ...) 
 res<-list(p.value=res$p.value, cliques=res$cliques)
 })
}
if (test=="var") {
xx<-lapply(pathways, function(p) {
cat(p[[2]],"\n")
 res<-CVT(x, group, p[[1]], alpha) 
 res<-list(p.value=res$p.value, cliques=res$cliques)
 })
 
 
 }
return(xx)
}
