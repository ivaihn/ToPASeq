rotate <-function (x, y, alpha, offset) 
{
    xn <- x * cos(alpha) - y * sin(alpha) + offset[1]
    yn <- x * sin(alpha) + y * cos(alpha) + offset[2]
    list(x = xn, y = yn)
}
drawHead<-function (type, xy, bbox, col, lwd, lty, len, out = TRUE) 
{
    db <- as.numeric(diff(bbox))
    dxy <- diff(xy) * db
    alpha <- atan(dxy[2]/dxy[1])
    r <- max(bbox)/130
    warn = FALSE
    normArrow <- function(r, alpha, xy, col, lwd, lty, out) {
        r <- r * 0.5
        x <- c(-1, 0, 1) * r
        y <- c(-1, 1, -1) * r
        off <- if (out) 
            90
        else -90
        alpha <- alpha - off * (pi/180)
        xyr <- rotate(x, y, alpha, xy[2, ])
        polygon(xyr, col = col, border = col, lwd = lwd, lty = lty)
    }
    switch(unlist(type), none = {
    }, box = {
        x <- c(-1, -1, 1, 1) * r
        y <- c(-1, 1, 1, -1) * r
        xyr <- rotate(x, y, alpha, xy[2, ])
        polygon(xyr, col = col, border = col, lwd = lwd, lty = lty)
    }, obox = {
        x <- c(-1, -1, 1, 1) * r
        y <- c(-1, 1, 1, -1) * r
        xyr <- rotate(x, y, alpha, xy[2, ])
        polygon(xyr, border = col, col = "white", lwd = lwd, 
            lty = lty)
    }, dot = {
        symbols(xy[2, 1], xy[2, 2], circles = r, inches = FALSE, 
            add = TRUE, fg = col, lwd = lwd, lty = lty, bg = col)
    }, odot = {
        symbols(xy[2, 1], xy[2, 2], circles = r, inches = FALSE, 
            add = TRUE, fg = col, lwd = lwd, lty = lty, bg = "white")
    }, diamond = {
        x <- c(-1, -1, 1, 1) * r
        y <- c(-1, 1, 1, -1) * r
        xyr <- rotate(x, y, alpha + 45 * (pi/180), xy[2, ])
        polygon(xyr, col = col, border = col, lwd = lwd, lty = lty)
    }, odiamond = {
        x <- c(-1, -1, 1, 1) * r
        y <- c(-1, 1, 1, -1) * r
        xyr <- rotate(x, y, alpha + 45 * (pi/180), xy[2, ])
        polygon(xyr, col = "white", border = col, lwd = lwd, 
            lty = lty)
    }, tee = {
        x <- c(0, 0) * r
        y <- c(-1, 1) * r
        xyr <- rotate(x, y, alpha, xy[2, ])
        lines(xyr, col = col, lwd = lwd, lty = lty)
    }, normal = {
        normArrow(r, alpha, xy, col, lwd, lty, out)
    }, open = {
        arrows(xy[1], xy[3], xy[2], xy[4], length = len, col = col, 
            lwd = lwd, lty = lty)
    }, vee = {
        arrows(xy[1], xy[3], xy[2], xy[4], length = len, col = col, 
            lwd = lwd, lty = lty)
    }, {
        warn <- TRUE
        arrows(xy[1], xy[3], xy[2], xy[4], length = len, col = col, 
            lwd = lwd, lty = lty)
    })
    warn
}
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

