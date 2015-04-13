
applyReduction<-function(reduct, NodeTable, agg.fun=mean) {
if (length(reduct)>0) {
r<-lapply(reduct, function(x) NodeTable[NodeTable[,3] %in% x,])
r<-t(sapply(r, function(x) sapply(1:10, function(i) 
 if (i==10) levels(x[,i])[max(as.numeric(x[,i]))] else 
 if (i==7) any(x[,i]) else
 if (i==5) agg.fun(x[,i]) else paste(x[,i], collapse=","))))
r<-data.frame(r)
names(r)<-names(NodeTable)
r$namesPlot<- names(reduct)

rownames(r)<-NULL
rbind(NodeTable,r)
} else NodeTable
}

edgetypeslegend<-function(x, cex){
par(mai=c(1.02,0,0,0))
plot(1, ylim=c(0,7), xlim=c(1,5), type="n", xlab="", ylab="", frame=FALSE, axes=FALSE)


text( 0.8, 6.5, "Edge types", pos=4, font=2, cex=cex)
lines(c(1,2), c(4,4), lty="dashed")
text( 2.2, 4, "dissociation", pos=4, cex=cex)
lines(c(1,2), c(4.5,4.5), lty="dashed")
lines(c(1.8,2,1.8), c(4.45,4.5,4.55), lty="solid")
text( 2.2, 4.5, "indirect effect", pos=4, cex=cex)

lines(c(1,2), c(5,5), lty="solid")
text( 2.2, 5, "bindig, association", pos=4, cex=cex)

lines(c(1,2), c(5.5,5.5), lty="solid")
lines(c(2,2), c(5.4,5.6), lty="solid")
text( 2.2, 5.5, "inhibition", pos=4, cex=cex)

lines(c(1,2), c(6,6), lty="solid")
lines(c(1.8,2,1.8), c(5.95,6,6.05), lty="solid")
text( 2.2, 6, "activation", pos=4, cex=cex)

text(0.8, 3.5, "Edge labels", pos=4, font=2, cex=cex)
text( 1.5, 3, "+p", cex=cex)
text( 2.2, 3, "phosphorylation", pos=4, cex=cex)
text( 1.5, 2.5, "-p", cex=cex)
text( 2.2, 2.5, "dephosphorylation", pos=4, cex=cex)

text( 1.5, 2, "+u", cex=cex)
text( 2.2, 2, "ubiquination", pos=4, cex=cex)

text(0.8, 1, "Edge colors", pos=4, font=2, cex=cex)

legend(x="bottom", legend=c("consistent", "inconstistent", "indecisive"), fill=c("black","red","grey87"), bty="n")
}

makeBreaks<-function(x, l, sym=TRUE){
if (sym) b<-c(seq(-max(abs(x)),0,length.out=l/2), seq(0,max(abs(x)),length.out=l/2)) else {}

b<-b[-which(b==0)[1]]
b
}

#vytvara tabulku s vlastnostami uzlov (mena (originalne, podla dat, pre kreslenie), pritomnost, genestats, farba vyplne, kategorie TSIG)
makeNodeTable<-function(g, gc, gp, breaks, deg.table, sigpal, tsig.whole, tsig, mis.col="white", p.th=0.05, col.lim=NULL ){
rownames(deg.table)<-as.character(deg.table[,1])
namesOrig<-nodes(g)
namesConv<-nodes(gc)
namesPlot<-nodes(gp)
isPresent<-nodes(gc) %in% as.character(deg.table[,1])
nodeStat <-deg.table[nodes(gc),2]
nodePval<-deg.table[nodes(gc),"pval"]
isSig<- nodePval<p.th

border<-ifelse(isSig | is.na(isSig),"white","black")

if (is.null(col.lim)) b<-makeBreaks(nodeStat[!is.na(nodeStat)] , breaks[1]) else
b<-makeBreaks(col.lim, breaks[1])
nodeCol  <-as.character(cut(nodeStat, b, labels=sigpal(length(b)-1), include.lowest=TRUE))
nodeCol[!isPresent]<-mis.col
#nodeCol[!isSig]<-mis.col

if (all(tsig.whole==1)) {tsigCat<-factor(tsig)} else {
b2<-makeBreaks(tsig.whole, breaks[2])
tsigCat<-cut(tsig, b2)
}


nodeTable<-data.frame(namesOrig=namesOrig, namesConv=namesConv, namesPlot=namesPlot,
                      isPresent=isPresent, nodeStat=nodeStat, nodePval=nodePval, isSig=isSig, border=border, nodeCol=nodeCol, tsigCat=tsigCat, stringsAsFactors=FALSE)
return(nodeTable)                      
}

#vytvara list s vlastnostami hran (sipky, ciarkovanie, znacky, farba)
makeEdgeList<-function(g, defaultEdgeAttrs){

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

undir<-c("binding","process","process(indirect)", "process(dissociation)")
tails<-setNames(rep("none", length(e.names)), e.names)
#tails[e.type %in% undir]<-"none"


dashed<-e.names[e.type %in% defaultEdgeAttrs[[1]][defaultEdgeAttrs[[1]][,2]  %in% c("indirect effect","dissociation"),1]]
dashed<-setNames(rep("dashed", length(dashed)), dashed)

edg.lab<-setNames(rep(NA, length(e.names)), e.names)
edg.lab[e.type %in% defaultEdgeAttrs[[1]][defaultEdgeAttrs[[1]][,2] == "phosphorylation",1]]<-"+p"
edg.lab[e.type %in% defaultEdgeAttrs[[1]][defaultEdgeAttrs[[1]][,2] == "dephosphorylation",1]]<-"-p"
edg.lab[e.type %in% defaultEdgeAttrs[[1]][defaultEdgeAttrs[[1]][,2] == "ubiquination",1]]<-"+u"
edg.lab<-edg.lab[!is.na(edg.lab)]

edg.col<-setNames(rep("black",nrow(edges(g))) , e.names )

EdgeList<-list(ENames=e.names, EType=e.type, Arrows=arrows, Tails=tails, Dashed=dashed, ELabels=edg.lab, Ecol=edg.col)
return(EdgeList)
}

#uprava hran a uzlov - farby podla toho ci platia, neplatia, ..
adjustAttr<-function(g, NodeTable, EdgeList,  stats, cols, remNodes=NULL){

if (is.character(NodeTable$isSig)) {
 sig<- NodeTable$isSig=="TRUE"
 names(sig)<-NodeTable$namesPlot
} else {
sig<-NodeTable$isSig
names(sig)<- nodes(g)
}
print(sig)
E<-edges(g)

selectValid <- selectInvalid<-rep(FALSE, nrow(E))

th<-0

selectValid[ stats[E[,1]] >  th & stats[E[,2]] >  th & (regexpr("activation",as.character(E[,4]))!=-1 | regexpr("expression",as.character(E[,4]))!=-1 )]<-TRUE
selectValid[ stats[E[,1]] >  th & stats[E[,2]] < -th & (regexpr("inhibition",as.character(E[,4]))!=-1 | regexpr("repression",as.character(E[,4]))!=-1)]<-TRUE
selectValid[ stats[E[,1]] < -th & stats[E[,2]] >  th & (regexpr("inhibition",as.character(E[,4]))!=-1 | regexpr("repression",as.character(E[,4]))!=-1)] <-TRUE
selectValid[ stats[E[,1]] < -th & stats[E[,2]] < -th & (regexpr("activation",as.character(E[,4]))!=-1 | regexpr("expression",as.character(E[,4]))!=-1)]<-TRUE
selectValid[!sig[E[,1] ] | !sig[E[,2] ] ]<-FALSE


selectInvalid[ stats[E[,1]] >  th & stats[E[,2]] >  th & (regexpr("inhibition",as.character(E[,4]))!=-1 | regexpr("repression",as.character(E[,4]))!=-1 )]<-TRUE
selectInvalid[ stats[E[,1]] >  th & stats[E[,2]] < -th & (regexpr("activation",as.character(E[,4]))!=-1 | regexpr("expression",as.character(E[,4]))!=-1)]<-TRUE
selectInvalid[ stats[E[,1]] < -th & stats[E[,2]] >  th & (regexpr("activation",as.character(E[,4]))!=-1 | regexpr("expression",as.character(E[,4]))!=-1)] <-TRUE
selectInvalid[ stats[E[,1]] < -th & stats[E[,2]] < -th & (regexpr("inhibition",as.character(E[,4]))!=-1 | regexpr("repression",as.character(E[,4]))!=-1)]<-TRUE
selectInvalid[!sig[E[,1] ] | !sig[E[,2] ] ]<-FALSE

valid.edges<-paste(E[,1], E[,2], sep="~")[selectValid]
not.valid.edges<-paste(E[,1], E[,2], sep="~")[selectInvalid]

valid.nodes<-unique(c(E[,1][selectValid | selectInvalid], E[,2][selectValid | selectInvalid]))
indecisive.edges<-paste(E[,1], E[,2], sep="~")[!selectValid & !selectInvalid] 

edg.col<-c(setNames(rep(cols[1],length(valid.edges)), valid.edges), 
           setNames(rep(cols[2],length(not.valid.edges)), not.valid.edges),
           setNames(rep(cols[3],length(indecisive.edges)), indecisive.edges) ) 
EdgeList$Ecol<-edg.col

if (!is.null(remNodes)) {
if (!is.character(remNodes)) stop("Invalid color specification for nodes with neutral interacions or with small expression change.")
valid.nodes<-valid.nodes[valid.nodes %in% NodeTable$namesPlot]
invalid.nodes<-as.character(NodeTable$namesPlot[ ! NodeTable$namesPlot %in% valid.nodes])
nodecol<-as.character(NodeTable$nodeCol)
nodecol[NodeTable$namesPlot %in% invalid.nodes]<-remNodes
NodeTable$nodeCol<-nodecol
}
out<-list(NodeTable, EdgeList)

return(out)
}

adjustAttrCli<-function(g, NodeTable, EdgeList, alpha, cliq, cols, th, remNodes=NULL){

E<-edges(g)
Edg<-paste(E[,1], E[,2],sep="~")
cliq<-cliq[alpha<th]

El<-lapply(cliq, function(x)
{
 if (length(x) > 2) 
  cli.edg<-c( apply(combn(x,2),2,          function(y) paste(y, collapse="~")), 
              apply(combn(x,2)[c(2,1),],2, function(y) paste(y, collapse="~"))) else
  if (length(x) == 2) cli.edg<-c(paste(x, collapse="~"), paste(x[c(2,1)], collapse="~")) else
  cli.edg<-character()
cli.edg<-cli.edg[cli.edg %in% Edg]
})
El<-El[sapply(El, length)>0]

setCols<-Map(function(x,y) {setNames(rep(y, length(x)),x)}, El, as.list(cols))

edg.col<-Reduce(c, setCols ) 

edg.col.other<-Edg[! Edg %in% names(edg.col)]
edg.col.other<-setNames(rep("grey87", length(edg.col.other)), edg.col.other)

EdgeList$Ecol<-c(edg.col, edg.col.other)


out<-list(NodeTable, EdgeList)

return(out)
}

renderOrig<-function(gp, NodeTable, EdgeList, nodesize, fontsize){
nA<-list(fillcolor = setNames(as.character(NodeTable$nodeCol), NodeTable$namesPlot),
         color = setNames(NodeTable$border, NodeTable$namesPlot), #"red" podla p-val
         height = setNames(nodesize + as.numeric(NodeTable$tsigCat), NodeTable$namesPlot),
         #fontsize = setNames(fontsize + as.numeric(NodeTable$tsigCat), NodeTable$namesPlot),
         label = setNames(as.character(NodeTable$namesPlot), NodeTable$namesPlot)
       
        )

eA<-list(
 label=setNames(EdgeList$ELabels,names(EdgeList$ELabels)),
 labelfontsize=setNames(rep(0.5*fontsize, length(EdgeList$ELabels)), names(EdgeList$ELabels)),
 color=setNames(EdgeList$Ecol,names(EdgeList$ELabels))
 )
xxg<-layoutGraph(as(gp,"graphNEL"), nodeAttrs=nA, edgeAttrs=eA)
nodeRenderInfo(xxg)<-list(fontsize = setNames(fontsize + as.numeric(NodeTable$tsigCat), NodeTable$namesPlot))
edgeRenderInfo(xxg)<-list(arrowhead=EdgeList$Arrows, arrowtail=EdgeList$Tails, lty=EdgeList$Dashed, col=EdgeList$Ecol)
return(xxg)
}

renderReduced<-function(gp, reduct, NodeTable, EdgeList, xxg, nodesize, fontsize, agg.fun){

nA<-list(fillcolor = setNames(as.character(NodeTable$nodeCol), NodeTable$namesPlot),
         color = setNames(NodeTable$border, NodeTable$namesPlot), #"red" podla p-val
         height = setNames(nodesize + as.numeric(NodeTable$tsigCat), NodeTable$namesPlot),
         #fontsize = setNames(fontsize + as.numeric(NodeTable$tsigCat), NodeTable$namesPlot),
         label = setNames(as.character(NodeTable$namesPlot), NodeTable$namesPlot)
        )

eA<-list(
 label=setNames(EdgeList$ELabels,names(EdgeList$ELabels)),
 labelfontsize=setNames(rep(0.5*fontsize, length(EdgeList$ELabels)), names(EdgeList$ELabels)),
 color=setNames(EdgeList$Ecol,names(EdgeList$ELabels))
 )

xxred<-layoutGraph(as(reduceGraph(gp, reduct),"graphNEL"), nodeAttrs=nA, edgeAttrs=eA )


if (length(reduct)>0) {
h<-sapply(reduct, function(x) ceiling(mean(nodeRenderInfo(xxg)$height[x])))
rad<-sapply(reduct, function(x) paste(nodeRenderInfo(xxg)$rWidth[x], collapse=","))
piecol<-sapply(reduct, function(x) paste(nA$fillcolor[x], collapse=","))
font<-sapply(reduct, function(x) ceiling(agg.fun(nodeRenderInfo(xxg)$fontsize[x])))

} else {h <- rad <- piecol <- font <- c()}


nodeRenderInfo(xxred)<-list(
 fontsize=c(setNames(fontsize + as.numeric(NodeTable$tsigCat), NodeTable$namesPlot),font), 
 height=c(nodeRenderInfo(xxg)$height,h), 
 radius=(c(nodeRenderInfo(xxg)$rWidth,  rad)),
 piecol=(c(nA$fillcolor[! names(nA$fillcolor) %in% unlist(reduct)], piecol ))
 )


edgeRenderInfo(xxred)<-list(arrowhead=EdgeList$Arrows, arrowtail=EdgeList$Tails, lty=EdgeList$Dashed, col=EdgeList$Ecol)
return(xxred)
}

getNodeCex<-function (g) 
{
    lw <- getRenderPar(g, "lWidth", "nodes")
    rw <- getRenderPar(g, "rWidth", "nodes")
    height <- getRenderPar(g, "height", "nodes")
    label <- unlist(getRenderPar(g, "label", "nodes"))
    if (is.null(label)) 
        label <- nodes(g)
 
    cex <- getRenderPar(g, "cex", "nodes")
    if (is.null(cex)) {
        nodeDims <- cbind(lw + rw, height)
        stw <- strwidth(label)
        sth <- strheight(label)
        strDims <- cbind(stw * 1.1, sth * 1.4)
        strDims[!nzchar(label), ] <- c(strwidth(" "), strheight(" "))
        cex <- min(nodeDims/strDims)
    }
    
   return(cex)
}

makeNodeSizeLegend<-function(bbox, categories, nodesize, fontsize, cex){


V<-categories
g.sub<-graphNEL(V,edgemode="directed")

if (length(V)==1) {
 nA.leg<-list(fontcolor=setNames("transparent","1"), 
              color=setNames("transparent","1"))} else 
 nA.leg<-list(height=c(nodesize+ setNames(seq_len(length(V)), V)),
              fontsize=c(fontsize + setNames(seq_len(length(V)), V))
 )

leg.xx<-layoutGraph(g.sub, nodeAttrs=nA.leg)
leg.xx@renderInfo@graph$bbox<-matrix(c( 0,0, round(bbox[2,1]*0.25,0),bbox[2,2]),2,2,  byrow=TRUE) 


n.leg<-length(nodes(leg.xx))
nodesx<-round(leg.xx@renderInfo@graph$bbox[2,1]/2,0) 
nodesy<-round(leg.xx@renderInfo@graph$bbox[2,2]/2/(n.leg+1)*(1:n.leg),0)
nodeswidth<-leg.xx@renderInfo@nodes$lWidth

leg.xx@renderInfo@nodes$rWidth<-nodeswidth
leg.xx@renderInfo@nodes$lWidth<-nodeswidth
leg.xx@renderInfo@nodes$height<-setNames(nodeswidth*2, nodes(leg.xx))
leg.xx@renderInfo@nodes$nodeX<-setNames(rep(nodesx, n.leg ), nodes(leg.xx))
leg.xx@renderInfo@nodes$nodeY<-setNames(nodesy, nodes(leg.xx))
leg.xx@renderInfo@nodes$labelX<-setNames(rep(nodesx+(2*max(nodeswidth)+0), n.leg ), nodes(leg.xx)) 
leg.xx@renderInfo@nodes$labelY<-setNames(nodesy, nodes(leg.xx))
leg.xx@renderInfo@nodes$cex<-setNames(rep(cex, n.leg ), nodes(leg.xx))

nodeRenderInfo(leg.xx)<-list(fontsize=c(fontsize + setNames(seq_len(length(V)), V)))


renderGraph(leg.xx, drawEdges=drawNoEdges<-function(g){} ,graph.pars=list(graph=list(main="")))
}

makeLegend<-function(bbox, categories, nodesize, fontsize, cex,  col.lim, br, sigpal, name, cex.legend){
makeNodeSizeLegend(bbox, categories, nodesize, fontsize, cex)

ext<-max(abs(col.lim))

fields::image.plot(add=TRUE, legend.only=TRUE, 
 zlim=c(-ext, ext),
 col=sigpal(br), legend.shrink=0.3,
 smallplot= c(.45,.55,0.6,0.9),
 legend.width=0.8, legend.lab=name, legend.mar=3.3, legend.line=3)

edgetypeslegend(1, cex.legend)
#ToPASeq:::edgetypeslegend(1)

}

headline<-function(res, which, plot){
if (is.character(which)) main<-which else 
if (is.list(res$topo.sig)) main<-names(res$topo.sig)[which] else
if (is.matrix(res$res)) main<-rownames(res$res)[which] else 
if (is.list(res$res)) main<-names(res$res)[which]


sig<-""
if (is.data.frame(res[[1]]) | is.matrix(res[[1]])) {
 if ("p" %in% colnames(res[[1]])) sig<- round(res$res[main,"p"],3)
 if ("p.value" %in% colnames(res[[1]])) sig<- round(res$res[main,"p.value"],3)
 if ("pG" %in% colnames(res[[1]])) sig<- round(res$res[main,"pG"],3)
 if ("alphaMean" %in% colnames(res[[1]])) sig<- round(res$res[main,"alphaMean"],3)
} 
if (is.list(res[[1]]) & length(res[[1]][[which]])==1) sig<-round(res$res[[which]][[1]]$p.value[2],3)

if (sig!="" & !is.na(sig)) sig<-paste("(p = ", sig,")", sep="") else sig=""
if (plot) title(paste(main, sig, sep="\n"), line=-2)
return(paste(main, sig, sep="\n"))
}

drawGraph<-function(xxred, res, which, NodeTable, nodesize, fontsize, statName,  cex.main=1, col.lim=NULL, breaks, legend=TRUE){


if (legend) layout(matrix(c(1,1,1,1,2,3), nrow=1)) 
 
renderGraph(xxred, graph.pars=list(graph=list(main=headline(res,which,FALSE), cex.main=cex.main)),  drawNodes=drawNodesPies2,drawEdges=renderEdgesTypes2)


bbox<-xxred@renderInfo@graph$bbox #4741 4561
categories<-levels(NodeTable$tsigCat)
CC<-getNodeCex(xxred)

if (is.null(col.lim)) {
geneNames<-as.character(NodeTable$namesConv)
isPr<-NodeTable$isPresent
isPr<-sapply(as.character(isPr), function(x) any(strsplit(x,",", fixed=TRUE)[[1]]=="TRUE"))
stats<-as.numeric(NodeTable$nodeStat)#res[[3]][geneNames[ isPr ]]
stats<-stats[!is.na(stats)]
col.lim<-range(stats)
}

if (legend) {makeLegend(bbox, categories, nodesize, fontsize, CC*1.1,  col.lim, breaks[1],"Log Fold-Change")
layout(1)
#par(.pardefault)
#graph.par(.graphpar)
}
}

#####