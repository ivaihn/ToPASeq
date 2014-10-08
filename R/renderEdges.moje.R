renderEdgesTypes<-function (g) 
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
        warn <- warn | suppressWarnings(renderSpline(splines[[i]], 
            arrowhead = arrowhead[i], arrowtail = arrowtail[i], 
            len = arrowLen, col = col[i], lty = lty[i], lwd = lwd[i], 
            bbox = getRenderPar(g, "bbox", "graph")))
    }
    if (warn) 
        warning("Unknown or unsupported arrowhead type. Using default instead.")
    TeachingDemos::shadowtext(labelX, labelY, label, col = textCol, cex = cex * as.numeric(fontsize)/14, 
 bg="white")
}

#getRenderPar<-Rgraphviz:::getRenderPar
#renderSpline<-Rgraphviz:::renderSpline

#renderGraph(xx, drawEdges=renderEdges.moje)

