pieRad<-function (x, xpos, ypos, labels = names(x), edges = 200, radius = 0.8, 
    density = NULL, angle = 45, col = NULL, border = NULL, lty = NULL, 
    main = NULL, ...) {
    if (!is.numeric(x) || any(is.na(x) | x <= 0)) 
        stop("pie: `x' values must be positive.")
    if (is.null(labels)) 
        labels <- as.character(1:length(x))
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    if (is.null(col)) 
        col <- if (is.null(density)) 
            c("lightblue", "mistyrose", "lightcyan", "lavender", 
                "cornsilk", "white")
        else par("fg")
    if (!is.null(col)) 
        col <- rep(col, length.out = nx)
    if (!is.null(border)) 
        border <- rep(border, length.out = nx)
    if (!is.null(lty)) 
        lty <- rep(lty, length.out = nx)
    if (!is.null(angle)) 
        angle <- rep(angle, length.out = nx)
    if (!is.null(density)) 
        density <- rep(density, length.out = nx)
    for (i in 1:nx) {
        n <- max(2, floor(edges * dx[i]))
        t2p <- 2 * pi * seq(x[i], x[i + 1], length = n)
        xc <- c(cos(t2p), 0) * radius[i] + xpos
        yc <- c(sin(t2p), 0) * radius[i] + ypos
        polygon(xc, yc, density = density[i], angle = angle[i], 
            border = NA, col = col[i], lty = lty[i])
       lines(xc[-length(xc)], yc[-length(yc)], col=border[i], lty = lty[i], lwd=1)
  
        
    }
    invisible(NULL)
}


drawNodesPies2<-function(g){
    nodeX <- getRenderPar(g, "nodeX", "nodes")
    nodeY <- getRenderPar(g, "nodeY", "nodes")
    lw <- getRenderPar(g, "lWidth", "nodes")
    rw <- getRenderPar(g, "rWidth", "nodes")
    height <- getRenderPar(g, "height", "nodes")
    labelX <- getRenderPar(g, "labelX", "nodes")
    labelY <- getRenderPar(g, "labelY", "nodes")
    radius <- getRenderPar(g, "radius", "nodes")
    piecol <- getRenderPar(g, "piecol", "nodes")
    fill <- unlist(getRenderPar(g, "fill", "nodes"))
    col <- unlist(getRenderPar(g, "col", "nodes"))
    lwd <- unlist(getRenderPar(g, "lwd", "nodes"))
    lty <- unlist(getRenderPar(g, "lty", "nodes"))
    textCol <- unlist(getRenderPar(g, "textCol", "nodes"))
    style <- unlist(getRenderPar(g, "style", "nodes"))
    shape <- getRenderPar(g, "shape", "nodes")
    label <- unlist(getRenderPar(g, "label", "nodes"))
    fontsize <- unlist(getRenderPar(g, "fontsize", "nodes"))
    if (is.null(label)) 
        label <- nodes(g)
    funs <- sapply(shape, is.function)
    possible.shapes <- c("circle", "ellipse", "box", "rectangle", 
        "plaintext", "triangle", "diamond")
    shape <- possible.shapes[pmatch(shape, possible.shapes, duplicates.ok = TRUE, 
        nomatch = 5)]
    i <- shape == "circle"
    if (any(i, na.rm = TRUE)) {
        rad <- pmin(height, (lw + rw))/2
        wh <- which(i)
        sapply(wh, function(ww) {
            r <- as.numeric(strsplit(as.character(radius[ww]), ",")[[1]])
            pc <- strsplit(as.character(piecol[ww]), ",")[[1]]
            b<- strsplit(as.character(col[ww]), ",")[[1]]
            if (any(b=="black")) 
             pieRad(rep(1, length(r)), nodeX[ww], nodeY[ww], edges = 200, radius = r, 
             angle = 45, col = pc, border = b, lty=2) else
             pieRad(rep(1, length(r)), nodeX[ww], nodeY[ww], edges = 200, radius = r, 
             angle = 45, col = pc, border = "transparent", lty=1)
        })
    }
    cex <- getRenderPar(g, "cex", "nodes")
    if (is.null(cex)) {
        nodeDims <- cbind(lw + rw, height)
        stw <- strwidth(label)
        sth <- strheight(label)
        strDims <- cbind(stw * 1.1, sth * 1.4)
        strDims[!nzchar(label), ] <- c(strwidth(" "), strheight(" "))
        cex <- min(nodeDims/strDims)
    }
    text(labelX, labelY, label, col = textCol, cex = cex * as.numeric(fontsize)/14)
   return(cex)
}

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

renderEdgesTypes2<-function (g) 
{
    lw <- getRenderPar(g, "lWidth", "nodes")
    rw <- getRenderPar(g, "rWidth", "nodes")
    height <- getRenderPar(g, "height", "nodes")
    splines <- getRenderPar(g, "splines", "edges")
    arrowhead <- unlist(getRenderPar(g, "arrowhead", "edges"))
    arrowtail <- getRenderPar(g, "arrowtail", "edges")
    label <- getRenderPar(g, "label", "edges")
    labelX <- getRenderPar(g, "labelX", "edges")
    labelY <- getRenderPar(g, "labelY", "edges")
    fontsize <- getRenderPar(g, "fontsize", "edges")
    textCol <- getRenderPar(g, "textCol", "edges")
    col <- unlist(getRenderPar(g, "col", "edges"))
    lty <- getRenderPar(g, "lty", "edges")
    lwd <- unlist(getRenderPar(g, "lwd", "edges"))
    cex <- getRenderPar(g, "cex", "edges")
    minDim <- min(max(rw + lw), max(height))
    arrowLen <- par("pin")[1]/diff(par("usr")[1:2]) * minDim/(1.5 * 
        pi)
    warn <- FALSE
    for (i in seq_along(splines)) {
       
if (length(splines[[i]])>1) {
n<-length(splines[[i]])
sp<-splines[[i]]
for ( j in seq_len(n))

if (j == 1)  {
suppressWarnings(renderSpline(sp[j], 
            arrowhead = "none", arrowtail = arrowtail[i], 
            len = arrowLen, col = col[i], lty = lty[i], lwd = lwd[i], 
            bbox = getRenderPar(g, "bbox", "graph")))
} else 
if (j == n)  { 
suppressWarnings(renderSpline(sp[j], 
            arrowhead = arrowhead[i], arrowtail = "none", 
            len = arrowLen, col = col[i], lty = lty[i], lwd = lwd[i], 
            bbox = getRenderPar(g, "bbox", "graph")))
} else { 
suppressWarnings(renderSpline(sp[j], 
            arrowhead = "none", arrowtail = "none", 
            len = arrowLen, col = col[i], lty = lty[i], lwd = lwd[i], 
            bbox = getRenderPar(g, "bbox", "graph")))
}
} else
suppressWarnings(renderSpline(splines[[i]], 
            arrowhead = arrowhead[i], arrowtail = arrowtail[i], 
            len = arrowLen, col = col[i], lty = lty[i], lwd = lwd[i], 
            bbox = getRenderPar(g, "bbox", "graph")))
}
    if (warn) 
        warning("Unknown or unsupported arrowhead type. Using default instead.")
    TeachingDemos::shadowtext(labelX, labelY, label, col = col, 
        cex = cex * as.numeric(fontsize)/14, bg = "white")
}
########