# x agopen vysledok drahy
# leg agopen vysledok uzlov pre legendu

plotRagraph<-function (x, y, leg, ...) {
    .local <- function (x, y, edgeAttrs = list(), ..., main = NULL, 
        cex.main = NULL, col.main = "black", sub = NULL, cex.sub = NULL, 
        col.sub = "black", drawNode = drawAgNode, xlab, ylab, 
        mai) 
    {
        if (missing(y)) 
            y <- x@layoutType
        x <- graphLayout(x, y)
        leg<-graphLayout(leg,y)
        plot.new()
        bg <- if (x@bg != "") 
            x@bg
        else par("bg")
        fg <- if (x@fg != "") 
            x@fg
        else par("fg")
        oldpars <- par(bg = bg, fg = fg)
        on.exit(par(oldpars), add = TRUE)
        if (missing(mai)) {
            mheight <- if (!is.null(main) && nchar(main) > 0) 
                sum(strheight(main, "inches", cex.main)) + 0.3
            else 0.1
            sheight <- if (!is.null(sub) && nchar(sub) > 0) 
                sum(strheight(main, "inches", cex.sub)) + 0.2
            else 0.1
            mai <- c(sheight, 0, mheight, 0)
        }
        oldpars <- par(mai = mai)
        on.exit(par(oldpars), add = TRUE)
        if (!is.null(sub) || !is.null(main)) 
            title(main, sub, cex.main = cex.main, col.main = col.main, 
                cex.sub = cex.sub, col.sub = col.sub)
        ur <- upRight(boundBox(x))
        bl <- botLeft(boundBox(x))
        ur.leg <- upRight(boundBox(leg))
        bl.leg <- botLeft(boundBox(leg))
       
        plot.window(xlim = c(getX(bl), getX(ur)), ylim = c(getY(bl), 
            getY(ur)), log = "", asp = NA, ...)
        usr<-par("usr")
        plot.window(xlim = c(getX(bl), getX(ur)*1.25), ylim = c(getY(bl), 
            getY(ur)), log = "", asp = NA, ...)
        coords<-list(xlim = c(getX(bl), getX(ur)*1.25), ylim = c(getY(bl), 
            getY(ur)))
        if (!missing(xlab) && !missing(ylab)) 
            stop("Arguments 'xlab' and 'ylab' are not handled.")
        agn <- AgNode(x)
        agn.leg<-AgNode(leg)
                
for (i in 1:length(agn.leg)) agn.leg[[i]]@center@x=getX(ur)*1.125#-2*getNodeRW(agn.leg[[length(agn.leg)]])
for (i in 1:length(agn.leg)) agn.leg[[i]]@center@y=i*(getY(ur)/(2*(length(agn.leg)+1)))
for (i in 1:length(agn.leg)) agn.leg[[i]]@txtLabel@labelColor="transparent" 

agn.leg2<-AgNode(leg)
                
for (i in 1:length(agn.leg2)) agn.leg2[[i]]@center@x=getX(ur)*1.125 + 2*getNodeRW(agn.leg[[length(agn.leg)]])
for (i in 1:length(agn.leg2)) agn.leg2[[i]]@center@y=i*(getY(ur)/(2*(length(agn.leg)+1)))
for (i in 1:length(agn.leg2)) agn.leg2[[i]]@color="transparent"

        nodeDims <- sapply(agn, function(n) {
            c(getNodeRW(n) + getNodeLW(n), getNodeHeight(n))
        })
        
        strDims <- sapply(agn, function(n) {
            s <- labelText(txtLabel(n))
            if (length(s) == 0) {
                rv <- c(strwidth(" "), strheight(" "))
            }
            else {
                rv <- c(strwidth(s) * 1.1, strheight(s) * 1.4)
            }
            return(rv)
        })
        cex <- min(nodeDims/strDims)
        if (is.finite(cex) && cex > 0) {
            oldpars <- par(cex = cex)
            on.exit(par(oldpars), add = TRUE)
        }
        if (length(drawNode) == 1) {
            lapply(agn, drawNode)
            lapply(agn.leg, drawNode)
            lapply(agn.leg2,drawNode)
        }
        else {
            if (length(drawNode) == length(AgNode(x))) {
                for (i in seq(along = drawNode)) {
                  drawNode[[i]](agn[[i]])
                }
            }
            else {
                stop(paste("Length of the drawNode parameter is ", 
                  length(drawNode), ", it must be either length 1 or the number of nodes.", 
                  sep = ""))
            }
        }
        arrowLen <- par("pin")[1]/(diff(usr[1:2])) * min(nodeDims)/pi
        lapply(AgEdge(x), lines, len = arrowLen, edgemode = edgemode, 
            ...)
        #invisible(x)
        invisible(coords)
    }
    .local(x, y, ...)
    

}

#coords<-plotRagraph(z,"dot", xx)