
text.lines<-function(x, which){
pos<-gregexpr(" ",x, fixed=TRUE)[[1]]
pos<-pos[seq(0, length(pos),which)]
letters<-strsplit(x,"")[[1]]
letters[pos]<-"\n"
x<-paste(letters,collapse="")
return(x)
}

.clique.edges<-function (nodes, pvalue) 
{
    edges <- expand.grid(nodes, nodes, stringsAsFactors = FALSE)
    diff <- edges[, 1] != edges[, 2]
    edges <- apply(edges[diff, ], 1, function(r) paste(r, collapse = "~"))
    pvalues <- rep(pvalue, length(edges))
    rbind(edges, pvalues)
}


plotCliques<-function(info, alpha = 0.05, color="red", node.color="white", nodesize, fontsize, add.legend=TRUE, layout="dot", intersp, ent=3) 
{
 
    if (length(edges(info$graph)) == 0) 
        stop("cannot render a graph with no edges")
    if (sum(!is.na(info$p.value)) == 0) {
        warning("No valid p-value, ploting general graph...")
        nnodes<-length(nodes(info$graph))
        g<-info$graph
        plot(g,  layout, nodeAttrs=list(fontsize=setNames(rep(fontsize,nnodes), nodes(g)), 
            width=setNames(rep(nodesize,nnodes), nodes(g))))
        }   else {
       
    g <- layoutGraph(info$graph)
    nnodes<-length(nodes(info$graph))
    check <- info$p.value < alpha
    check <- check[!is.na(check)]
    significant <- info$cliques[check]
    pvalues <- info$p.value[check]
    if (length(significant)) {
        nodes <- unique(unlist(significant))
        edges.rd <- matrix(data = unlist(sapply(1:length(significant), 
            function(i) .clique.edges(significant[[i]], pvalues[[i]]), 
            simplify = FALSE)), ncol = 2, byrow = TRUE)
        edges <- as.matrix(tapply(edges.rd[, 2], edges.rd[, 1], 
            function(x) min(as.numeric(x))))
             
      if (is.character(color)) palette<-rep(color, sum(check)) else
      if (is.function(color)) palette<-color(sum(check)) else
       stop("'color' must be either character or function")
      
            colors<- palette[as.numeric(factor(edges))]
            colors[is.na(colors)]<-"grey"
            leg.names<-sapply(info$cliques[match(levels(factor(edges)), info$p.value)], function(x) paste(x, collapse=" "))
            leg.names<-paste("Nodes: ", leg.names, " (p=",round(info$p.value[match(levels(factor(edges)), info$p.value)], 3), ")", sep="")
            leg.names<-unlist(lapply(leg.names, function(x) text.lines(x, ent)))
            names(colors) <- rownames(edges)
            
        
        edgeAttrs<- list(color = colors, weight=setNames(rep(2, nrow(edges)), rownames(edges)))
        colors <- rep(node.color, length(nodes))
        names(colors) <- nodes
        nodeAttrs<- list(fillcolor = colors)
    }
    
    if (add.legend) {
    layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
    plot(g,layout, edgeAttrs=edgeAttrs, nodeAttrs=c(nodeAttrs, 
    list(fontsize=setNames(rep(fontsize,nnodes), nodes(g)), 
            width=setNames(rep(nodesize,nnodes), nodes(g)))))
    
    if (length(significant)) {
       
            par(mai=c(0,0,0,0))
            plot(1,1,  type="n", xlab="", ylab="", axes=FALSE, frame.plot=FALSE, main="");
            legend(x="center", bty="n", pt.cex=2,pch=15, legend=leg.names, col=palette, y.intersp=intersp )
     }            
     layout(1)} else
    plot(g,layout, edgeAttrs=edgeAttrs, nodeAttrs=c(nodeAttrs, 
    list(fontsize=setNames(rep(fontsize,nnodes), nodes(g)), 
            width=setNames(rep(nodesize,nnodes), nodes(g)))))
}    

 }

########################
plot.topResult<-function(x, which, graphs, stats="logFC", convert=TRUE, IDs="entrez", graphIDs="symbol", col.lim=NULL,reduction=list(), agg.fun=function(x) mean(x, na.rm=TRUE),
 logical=NULL, sig.th=0.1, title=TRUE, cex.main, breaks=c(100,5),
  pallete.colors=c("blue","white", "red"), na.col="grey87", cli.color="red", layout="dot", nodesize=1, fontsize=14, 
  alpha=0.05, add.legend=TRUE, statName="Log fold change", ...  ){



#if (!require("Rgraphviz")) stop("Rgraphviz package is missing, please install it")
res<-x

 g<-graphs[[which]]
 if (convert) gc<-convertIdentifiers(graphs[[which]],IDs) else gc<-g
 gp<-convertIdentifiers(graphs[[which]],graphIDs)
 deg.table<-x$degtable
 
 sigpal<-colorRampPalette(pallete.colors)
 na.col<-colorRampPalette(c(na.col))(1)

 defaultEdgeAttrs<-makeDefaultEdgeData()
 
if ("topResultC" %in% class(res)) {
 cliq<-res$topo.sig[[which]]
 NodeTable<-makeNodeTable(g, gc, gp, breaks, deg.table, sigpal, tsig.whole, tsig, mis.col=na.col, p.th=alpha, col.lim=col.lim )
 EdgeList<-makeEdgeList(gp, defaultEdgeAttrs)

cols<-cli.color
 att<-adjustAttrCli(gc, NodeTable, EdgeList, cliq[[1]], cliq[[2]], cols, alpha, remNodes=NULL)

xxg<-renderOrig(gp, NodeTable, EdgeList, nodesize, fontsize)
xxred<-renderReduced( gp, reduction, att[[1]], att[[2]], xxg, fontsize)
drawGraph(xxred, res, which, NodeTable, nodesize, fontsize, statName=statName, cex.main=cex.main, legend=add.legend)
}

if ("topResultW" %in% class(res) ){
tsig<-res$topo.sig[[which]][1,]; 
 if (length(tsig)==0) stop("This pathway was not analysed") 
 tsig.whole<-unlist(sapply(res$topo.sig, function(x) x[1,]))
 }
if ("topResultE" %in% class(res) ) {
nod<-nodes(graphs[[which]])
tsig<-setNames(rep(1, length(nod)), nod)
tsig.whole<-rep(1, 100)
}
 
if ("topResultW" %in% class(res) | "topResultE" %in% class(res) ) {
 
NodeTable<-makeNodeTable(g, gc, gp, breaks,  deg.table, sigpal, tsig.whole, tsig, col.lim=col.lim)
 EdgeList<-makeEdgeList(gp, defaultEdgeAttrs)
 if (length(reduction)>0) {gpr<-reduceGraph(gp, reduction)} else gpr<-gp
 
 NodeTable.red<-applyReduction(reduction, NodeTable, agg.fun)
 EdgeList.red<-makeEdgeList(gpr, defaultEdgeAttrs)
 
 if (logical) {
 stats<-setNames(as.numeric(NodeTable.red$nodeStat), NodeTable.red$namesPlot)
 att<-adjustAttr(gpr, NodeTable.red, EdgeList.red, stats, cols=c("black","red", "grey87"), remNodes="white")
 } else att<-list(NodeTable, EdgeList)
 
 xxg<-renderOrig(gp, NodeTable, EdgeList, nodesize, fontsize)
 xxred<-renderReduced( gp, reduction, att[[1]], att[[2]], xxg, fontsize)
 drawGraph(xxred, res, which, NodeTable, nodesize, fontsize, statName=statName, cex.main=cex.main, col.lim=col.lim, legend=add.legend)
}




}