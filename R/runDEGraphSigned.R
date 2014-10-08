runDEGraphSigned<-function (pathways, expr, classes, useInteractionSigns=FALSE, both.directions, maxNodes = 150)
{
    if (!require(DEGraph))
        stop("library DEGraph is missing")

    checkPkgVersion<-function (name, min_version)
{
    version <- package_version(installed.packages()[name, "Version"])
    if (version < package_version(min_version))
        stop("the installed ", name, " version is too old (need at least ",
            min_version, ")")
}

    checkPkgVersion("DEGraph", "1.4.0")


    pathways <- filterPathwaysByNodeNum(pathways, maxNodes)
    lapplyCapturingErrors(pathways, function(p) {
    if (insufficientCommonGenes(p, rownames(expr)))
        return(NULL)
    if (both.directions) g <- buildGraphNEL(nodes(p), edges(p), TRUE) else
          g <- buildGraphNEL(nodes(p), edges(p), FALSE)
    #g <-  pathwayGraph(p)
    if (useInteractionSigns) {

    eA<-merge(defaultEdgeAttrs[[1]], defaultEdgeAttrs[[2]], by.x=2, by.y=1, all=TRUE)
    pos<-as.character(eA[,2][!is.na(eA[,2]) & eA[,"beta"]==1])
    neg<-as.character(eA[,2][!is.na(eA[,2]) & eA[,"beta"]==-1])
    neu<-as.character(eA[,2][!is.na(eA[,2]) & eA[,"beta"]==0])


    signMat<-as(g,"matrix")
    signMat[,]<-0

    e<-edges(p)
    posedg<-e[,1:2][e[,4] %in% pos, ]
    posind<-nrow(signMat)*(match(posedg[,2], colnames(signMat))-1)+match(posedg[,1], colnames(signMat))
    signMat[posind] <-  1
    negedg<-e[,1:2][e[,4] %in% neg, ]
    negind<-nrow(signMat)*(match(negedg[,2], colnames(signMat))-1)+match(negedg[,1], colnames(signMat))
    signMat[negind] <- -1

    neuedg<-e[,1:2][e[,4] %in% neu, ]
    g<-removeEdge(from=neuedg[,1], to=neuedg[,2], g)
    g@graphData$signMat<-signMat
    #print(g)

    testOneGraph(g, expr, classes, useInteractionSigns = FALSE)

    } else {
    #print(g)

    testOneGraph(g, expr, classes, useInteractionSigns = FALSE)
    }

})
}