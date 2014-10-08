 drawNodesPies<-function(g){

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
             r<-as.numeric(strsplit(as.character(radius[ww]),",")[[1]])
             pc<-strsplit(as.character(piecol[ww]),",")[[1]]
             pieRad(rep(1,length(r)), nodeX[ww], nodeY[ww], edges = 200, radius = r, angle = 45, col = pc, border="transparent")

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

}

#set.seed(123)
# V <- letters[1:10]
# M <- 1:4
# g1 <- randomGraph(V, M, 0.2)
# g1layout <- layoutGraph(g1)
#renderGraph(g1layout, drawNodes=drawNodesPies, main="Example Pie Chart Plot")

