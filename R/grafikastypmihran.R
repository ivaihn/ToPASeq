plotGraphW<-function(res, which, graphs, convert, IDs, graphIDs, reduction, combnfunction, logical, sig.th, pallete.colors, na.col,breaks, layout, nodesize, fontsize, shift, tsig, tsig.whole, stats, title){
.pardefault <- par(no.readonly = TRUE)
.graphpar<-graph.par()

if (convert) g<-convertIdentifiers(graphs[[which]],IDs) else
   g<-graphs[[which]]

g.orig<-g
pres.nodes<-nodes(g)[nodes(g) %in% names(shift)]
miss.nodes<-nodes(g)[!nodes(g) %in% names(shift)]


if (length(reduction)!=0) {
g<-reduceGraph(g, reduction)
tsig<-c(tsig, sapply(reduction, function(x, res, combnfunction) combnfunction(res[[2]][[which]][1,x]), res, combnfunction))
shift<-c(shift, sapply(reduction, function(x, res, combnfunction) combnfunction(res[[3]][x]), res, combnfunction))
}

gr<-pathwayGraph(g)
#palletes
sigpal<-colorRampPalette(pallete.colors)
na.col<-colorRampPalette(c(na.col))(1)

# color matrix


b<-c(seq(-max(abs(shift)),0,length.out=breaks[1]/2), seq(0,max(abs(shift)),length.out=breaks[1]/2))
b<-b[-which(b==0)[1]]
b.shift<-b
colp1<-cut(shift, b, labels=sigpal(length(b)-1), include.lowest=TRUE)
names(colp1)<-names(shift)


if (all(tsig.whole==1)) {tsig.cut<-factor(tsig)} else { 
b<-hist(tsig.whole, breaks=breaks[2], plot=FALSE)$breaks
tsig.cut<-cut(tsig, b, include.lowest=TRUE)
names(tsig.cut)<-names(tsig) 
}

#graph<-convertIdentifiers(graphs[[which]], graphIDs)
if (length(reduction)==0 & convert) node.names<-nodes(convertIdentifiers(graphs[[which]], graphIDs))  else   
  node.names<-nodes(g)
   
nA<-list(fillcolor=setNames(as.character(colp1[pres.nodes]), pres.nodes), 
         color=setNames(rep("white", length(nodes(gr))), nodes(gr)),
         height=nodesize+c(setNames(as.numeric(tsig.cut[pres.nodes]), pres.nodes), setNames(rep(0.8,length(miss.nodes) ), miss.nodes)),
         fontsize=fontsize+c(setNames(as.numeric(tsig.cut[pres.nodes]), pres.nodes), setNames(rep(0.8,length(miss.nodes) ), miss.nodes)),
         label=setNames(node.names, nodes(gr))  
    )


#setting up edge attributes
e.names<-edges(g)[,1:2]
e.names<-paste(e.names[,1], e.names[,2], sep="~")
e.type<-as.character(edges(g)[,4])

openarr<-c("activation", "indirect effect", "dephosphorylation", "phosphorylation", "ubiquination")
noarr<-c("binding/association", "dissociation")
arrows<-setNames(rep(NA, length(e.names)), e.names)
arrows[e.type %in% defaultEdgeAttrs[[1]][defaultEdgeAttrs[[1]][,2] %in% openarr,1]]<-"open"
arrows[e.type %in% defaultEdgeAttrs[[1]][defaultEdgeAttrs[[1]][,2] %in% noarr,1]]<-"none"
arrows[e.type %in% defaultEdgeAttrs[[1]][defaultEdgeAttrs[[1]][,2] == "inhibition",1]]<-"tee"
arrows<-arrows[!is.na(arrows)]



tails<-setNames(rep(NA, length(e.names)), e.names)
#tails[e.type %in% defaultEdgeAttrs[[1]][defaultEdgeAttrs[[1]][,2] %in% noarr,1]]<-"none"
tails<-tails[!is.na(tails)]



dashed<-e.names[e.type %in% defaultEdgeAttrs[[1]][defaultEdgeAttrs[[1]][,2]  %in% c("indirect effect","dissociation"),1]]
dashed<-setNames(rep("dashed", length(dashed)), dashed)

edg.lab<-setNames(rep(NA, length(e.names)), e.names)
edg.lab[e.type %in% defaultEdgeAttrs[[1]][defaultEdgeAttrs[[1]][,2] == "phosphorylation",1]]<-"+p"
edg.lab[e.type %in% defaultEdgeAttrs[[1]][defaultEdgeAttrs[[1]][,2] == "dephosphorylation",1]]<-"-p"
edg.lab[e.type %in% defaultEdgeAttrs[[1]][defaultEdgeAttrs[[1]][,2] == "ubiquination",1]]<-"+u"
edg.lab<-edg.lab[!is.na(edg.lab)]


if (is.null(logical)) edg.col<-setNames(rep("black",nrow(edges(g))) , e.names ) else 
if (logical) {               
th=sig.th
log.edg<-edges(g)
tmp<-rep(FALSE, nrow(log.edg))
tmp[ res[[3]][log.edg[,1]] >  th & res[[3]][log.edg[,2]] >  th & (regexpr("activation",as.character(log.edg[,4]))!=-1 | regexpr("expression",as.character(log.edg[,4]))!=-1 )]<-TRUE
tmp[ res[[3]][log.edg[,1]] >  th & res[[3]][log.edg[,2]] < -th & (regexpr("inhibition",as.character(log.edg[,4]))!=-1 | regexpr("repression",as.character(log.edg[,4]))!=-1)]<-TRUE
tmp[ res[[3]][log.edg[,1]] < -th & res[[3]][log.edg[,2]] >  th & (regexpr("inhibition",as.character(log.edg[,4]))!=-1 | regexpr("repression",as.character(log.edg[,4]))!=-1)] <-TRUE
tmp[ res[[3]][log.edg[,1]] < -th & res[[3]][log.edg[,2]] < -th & (regexpr("activation",as.character(log.edg[,4]))!=-1 | regexpr("expression",as.character(log.edg[,4]))!=-1)]<-TRUE

#all.edg.names<-paste(log.edg[,1], log.edg[,2], sep="~")
#ok.edg.names<-paste(log.edg[,1], log.edg[,2], sep="~")[tmp]
#not.ok.edg<-setdiff(all.edg.names, ok.edg.names)
not.ok.edg<-paste(log.edg[,1], log.edg[,2], sep="~")[!tmp]
edg.col<-setNames(rep("gray87",length(not.ok.edg)), not.ok.edg)

nA$fillcolor[!( names(nA$fillcolor) %in% unique(c(log.edg[,1][tmp], log.edg[,2][tmp])))]<-"white"
 } else   {
th=sig.th
log.edg<-edges(g)
tmp<-rep(FALSE, nrow(log.edg))
tmp[ res[[3]][log.edg[,1]] >  th & res[[3]][log.edg[,2]] >  th & (regexpr("inhibition",as.character(log.edg[,4]))!=-1 | regexpr("repression",as.character(log.edg[,4]))!=-1 )]<-TRUE
tmp[ res[[3]][log.edg[,1]] >  th & res[[3]][log.edg[,2]] < -th & (regexpr("activation",as.character(log.edg[,4]))!=-1 | regexpr("expression",as.character(log.edg[,4]))!=-1)]<-TRUE
tmp[ res[[3]][log.edg[,1]] < -th & res[[3]][log.edg[,2]] >  th & (regexpr("activation",as.character(log.edg[,4]))!=-1 | regexpr("expression",as.character(log.edg[,4]))!=-1)] <-TRUE
tmp[ res[[3]][log.edg[,1]] < -th & res[[3]][log.edg[,2]] < -th & (regexpr("inhibition",as.character(log.edg[,4]))!=-1 | regexpr("repression",as.character(log.edg[,4]))!=-1)]<-TRUE

#all.edg.names<-paste(log.edg[,1], log.edg[,2], sep="~")
#ok.edg.names<-paste(log.edg[,1], log.edg[,2], sep="~")[tmp]
#not.ok.edg<-setdiff(all.edg.names, ok.edg.names)
not.ok.edg<-paste(log.edg[,1], log.edg[,2], sep="~")[!tmp]
edg.col<-setNames(rep("gray87",length(not.ok.edg)), not.ok.edg)

nA$fillcolor[!( names(nA$fillcolor) %in% unique(c(log.edg[,1][tmp], log.edg[,2][tmp])))]<-"white" 
 }        
 

xxg<-layoutGraph(pathwayGraph(g.orig), nodeAttrs=nA[c("color","height","label","fontsize")], edgeAttrs=list(label=edg.lab, labelfontsize=setNames(rep(fontsize+10, length(edg.lab)), names(edg.lab)) ))
nodeRenderInfo(xxg)<-list(fontsize=c(fontsize+c(setNames(as.numeric(tsig.cut[pres.nodes]), pres.nodes), setNames(rep(0.8,length(miss.nodes) ), miss.nodes))))
              

if (length(reduction)>0){
h<-sapply(reduction, function(x) ceiling(combnfunction(nA$height[x])))
nA$height=c(nA$height, h)
font<-sapply(reduction, function(x) ceiling(combnfunction(nodeRenderInfo(xxg)$fontsize[x])))
nA$fontsize=c(nA$fontsize, font)
} else {h<-c(); font<-c()}


xx<-layoutGraph(gr, nodeAttrs=nA[c("color","height", "label","fontsize")], edgeAttrs=list(label=edg.lab, labelfontsize=setNames(rep(fontsize+10, length(edg.lab)), names(edg.lab)) ))
edgeRenderInfo(xx)<-list(arrowhead=arrows, arrowtail=tails, lty=dashed, col=edg.col) #, textCol=setNames(rep("red", length(edg.lab)), names(edg.lab)) , labelX=edgeRenderInfo(xx)$labelX+20)



#rw<-sapply(reduction, function(x) ceiling(combnfunction(nodeRenderInfo(xxg)$rWidth[x])))
#lw<-sapply(reduction, function(x) ceiling(combnfunction(nodeRenderInfo(xxg)$lWidth[x])))
if (length(reduction)>0) {
h<-sapply(reduction, function(x) ceiling(combnfunction(nodeRenderInfo(xxg)$height[x])))

rad<-sapply(reduction, function(x) paste(nodeRenderInfo(xxg)$rWidth[x], collapse=","))
piecol<-sapply(reduction, function(x) paste(nA$fillcolor[x], collapse=","))
} else {h<-c(); rad<-c(); piecol<-c()}


nodeRenderInfo(xx)<-list(fontsize=c(fontsize+c(setNames(as.numeric(tsig.cut[pres.nodes]), pres.nodes), setNames(rep(0.8,length(miss.nodes) ), miss.nodes)),font), 
                         height=c(nodeRenderInfo(xx)$height,h), #lWidth=lw, rWidth=rw,
                         radius=(c(setNames(nodeRenderInfo(xx)$rWidth[setdiff(pres.nodes, unlist(reduction))], setdiff(pres.nodes, unlist(reduction))),  rad)),
                         piecol=(c(nA$fillcolor[! names(nA$fillcolor) %in% unlist(reduction)], piecol ))
                         #height=nodesize+c(setNames(as.numeric(tsig.cut[pres.nodes]), pres.nodes), setNames(rep(0.8,length(miss.nodes) ), miss.nodes))
                        )
             
#if (title) {headline(res,which)}
#par.def<-graph.par()

layout(matrix(c(1,1,1,1,2,3), nrow=1)) 

#xx@renderInfo@graph$bbox[2,1]<- xx@renderInfo@graph$bbox[2,1]*1.25
#suppressWarnings()
renderGraph(xx, graph.pars=list(graph=list(main=headline(res,which,FALSE), cex.main=1.8)), drawNodes=drawNodesPies2, drawEdges=renderEdgesTypes2)


bbox<-xx@renderInfo@graph$bbox


V<-levels(tsig.cut)
g.sub<-graphNEL(V,edgemode="directed")

if (length(V)==1) {
 nA.leg<-list(fontcolor=setNames("transparent","1"), 
              color=setNames("transparent","1"))} else 
 nA.leg<-list(height=c(nodesize+ setNames(as.numeric(1:nlevels(tsig.cut)), V)),
              fontsize=c(fontsize + setNames(as.numeric(1:nlevels(tsig.cut)), V))
 )

leg.xx<-layoutGraph(g.sub, nodeAttrs=nA.leg)

leg.xx@renderInfo@graph$bbox<-matrix(c( 0,0, round(bbox[2,1]*0.25,0),bbox[2,2]),2,2,  byrow=TRUE) 


n.leg<-length(nodes(leg.xx))
nodesx<-round(leg.xx@renderInfo@graph$bbox[2,1]/2,0) 
nodeswidth<-leg.xx@renderInfo@nodes$lWidth

leg.xx@renderInfo@nodes$rWidth<-nodeswidth
leg.xx@renderInfo@nodes$lWidth<-nodeswidth
leg.xx@renderInfo@nodes$height<-setNames(nodeswidth*2, nodes(leg.xx))
leg.xx@renderInfo@nodes$nodeX<-setNames(rep(nodesx, n.leg ), nodes(leg.xx))
leg.xx@renderInfo@nodes$nodeY<-setNames(round(leg.xx@renderInfo@graph$bbox[2,2]/2/(n.leg+1)*(1:n.leg),0), nodes(leg.xx))
leg.xx@renderInfo@nodes$labelX<-setNames(rep(nodesx+(1.5*max(nodeswidth)+0), n.leg ), nodes(leg.xx)) #max(2*leg.xx@renderInfo@nodes$labelWidth)+2
leg.xx@renderInfo@nodes$labelY<-setNames(round(leg.xx@renderInfo@graph$bbox[2,2]/2/(n.leg+1)*(1:n.leg),0), nodes(leg.xx))
#nodeRenderInfo(leg.xx)<-list(fontsize=mean(nodeRenderInfo(xx)$fontsize))
nodeRenderInfo(leg.xx)<-list(fontsize=c(fontsize + setNames(as.numeric(1:nlevels(tsig.cut)), V)))

#nodeRenderInfo(leg.xx)<-list(fontsize=fontsize)

renderGraph(leg.xx, drawEdges=drawNoEdes<-function(g){} ,graph.pars=list(graph=list(main="")))


fields::image.plot(add=TRUE, legend.only=TRUE, 
 zlim=range(b.shift), 
 col=sigpal(length(b.shift)-1), legend.shrink=0.3,
 smallplot= c(.45,.55,0.6,0.9),
 legend.width=0.8, legend.lab=stats, legend.mar=3.3, legend.line=3)


edgetypeslegend(1)
par(.pardefault)
graph.par(.graphpar)
}
########################
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