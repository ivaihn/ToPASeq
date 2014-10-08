convertIDs<-function(pathway, ID="entrez"){
out<-sapply(pathway, convertIdentifiers, ID)
return(out)
}

#SPIA
SPIA<-function(x, group, pathways, type="MA", convert=TRUE,  IDs="entrez",  gene.stat="logFC",both.directions=TRUE,logFC.th=2, p.val.th=0.05,test=NULL, edgeAttrs=defaultEdgeAttrs,  ...)
{
if (any(is(x) == "ExpressionSet")) {
exprs<-exprs(x) 
if (is.character(group) & ! any(group == colnames(pData(x)))) stop( paste( group,"is not present in phenoData"))  else
if (is.numeric(group)  & group > ncol(pData(x))) stop(paste("phenoData contain only", ncol(pData(x)), "columns")) else
if (is.numeric(group) | (is.character(group)) ) group<-pData(x)[,group] else 
 stop ("'group' is must be either numeric or character")
 }  else exprs<-x

if (!is.list(pathways) | is(pathways[[1]])!="pathway") stop("Argument 'pathways' is not a list of pathway-class objects")
edgelevels<-levels(unlist(sapply(pathways, function(x) edges(x)[,4])))

if (all(edgelevels=="process(indirect effect")) warning("All interactions are set to 'process(indirect effect)'. Pathways probably do not contain valid information about interation type.")
if (convert) {CheckNames(convertIDs(pathways,IDs), exprs)} else CheckNames(pathways,exprs)

if (type=="MA") out<-runSPIA(pathways, exprs, group,  gene.stat, both.directions, convert=convert, logFC.th=logFC.th, p.val.th=p.val.th, IDs=IDs, edgeAttrs=edgeAttrs, ...)
if (type=="RNASeq") {
if (is.null(test)) {
  test<-"vstlimma"
  cat("test was not specified. 'vstlimma' used as default\n")
  out<-runSPIA4RNASeq(pathways, exprs, group, gene.stat, both.directions=both.directions, convert=convert, logFC.th=logFC.th, p.val.th=p.val.th, IDs=IDs, test=test, edgeAttrs=edgeAttrs,...)
} else 
out<-runSPIA4RNASeq(pathways, exprs, group, gene.stat, both.directions=both.directions, convert=convert, logFC.th=logFC.th, p.val.th=p.val.th, IDs=IDs, test=test, edgeAttrs=edgeAttrs, ...)
}
return(out)
}

#TAPPA
TAPPA<-function(x, group, pathways, type="MA", convert=TRUE, IDs="entrez", gene.stat="logFC", both.directions=TRUE, normalize=TRUE, verbose=FALSE, norm.method=NULL)
{

if (any(is(x) == "ExpressionSet")) {
exprs<-exprs(x) 
if (is.character(group) & ! any(group == colnames(pData(x)))) stop( paste( group,"is not present in phenoData"))  else
if (is.numeric(group)  & group > ncol(pData(x))) stop(paste("phenoData contain only", ncol(pData(x)), "columns")) else
if (is.numeric(group) | (is.character(group)) ) group<-pData(x)[,group] else 
 stop ("'group' is must be either numeric or character")
 }  else exprs<-x

if (!is.list(pathways) | is(pathways[[1]])!="pathway") stop("Argument 'pathways' is not a list of pathway-class objects")
if (convert) pathways<-convertIDs(pathways, IDs)
CheckNames(pathways, exprs)

if (type=="MA") out<-runTAPPA(pathways, exprs, group,gene.stat, both.directions, normalize, verbose)
if (type=="RNASeq") {
 if (is.null(norm.method)) {
  norm.method<-"TMM"
  cat("Normalization method was not specified. TMM used as default\n")
  out<-runTAPPA4RNASeq(pathways, exprs, group,gene.stat, both.directions, normalize, verbose, norm.method)
  } else
  out<-runTAPPA4RNASeq(pathways, exprs, group,gene.stat, both.directions, normalize, verbose, norm.method)
  }

return(out)
}
#DEGraph
DEGraph<-function(x, group, pathways, type="MA", convert=TRUE, IDs="entrez",gene.stat="logFC",both.directions=TRUE, overall="biggest",  useInteractionSigns=TRUE, norm.method=NULL)
{
if (any(is(x) == "ExpressionSet")) {
expr<-exprs(x)
if (is.character(group) & ! any(group == colnames(pData(x)))) stop( paste( group,"is not present in phenoData"))  else
if (is.numeric(group)  & group > ncol(pData(x))) stop(paste("phenoData contain only", ncol(pData(x)), "columns")) else
if (is.numeric(group) | (is.character(group)) ) group<-pData(x)[,group] else 
 stop ("'group' is must be either numeric or character")
 }  else expr<-x
  
  
if (!is.list(pathways) | is(pathways[[1]])!="pathway") stop("Argument 'pathways' is not a list of pathway-class objects")
if (convert) pathways<-convertIDs(pathways, IDs)
CheckNames(pathways, expr)

if (type=="MA") {
  res<-runDEGraphSigned(pathways, expr, group, useInteractionSigns, both.directions, maxNodes = 150)
  res<-processDEGraph(res[[1]],overall)
  deg.table<-testlimma(expr, group)
  if (gene.stat=="logFC") degt<-deg.table$logFC
  if (gene.stat=="mod.T") degt<-deg.table$t
  names(degt)<-deg.table$ID
  out<-list(res=res, topo.sig=NULL, degtest=degt)
  }

if (type=="RNASeq") {
 if (is.null(norm.method)) {
  norm.method<-"TMM"
  cat("Normalization method was not specified. TMM used as default\n")
  out<-runDEGraph4RNASeq(pathways, expr, group,gene.stat, overall,  norm.method=norm.method, useInteractionSigns=useInteractionSigns, both.directions=both.directions)
  } else
  out<-runDEGraph4RNASeq(pathways, expr, group,gene.stat, overall,  norm.method=norm.method,useInteractionSigns=useInteractionSigns, both.directions=both.directions)
  }
class(out)<-c(class(out), "topResultE","topResult")
return(out)
}
#TopologyGSA
TopologyGSA<-function(x, group, pathways, type="MA", convert=TRUE, IDs="entrez", both.directions=TRUE, test="mean", testCliques=FALSE, alpha=0.05, nperm=10000, norm.method=NULL, maxNodes=150 )
{
if (any(is(x) == "ExpressionSet")) {
exprs<-exprs(x) 
if (is.character(group) & ! any(group == colnames(pData(x)))) stop( paste( group,"is not present in phenoData"))  else
if (is.numeric(group)  & group > ncol(pData(x))) stop(paste("phenoData contain only", ncol(pData(x)), "columns")) else
if (is.numeric(group) | (is.character(group)) ) group<-pData(x)[,group] else 
 stop ("'group' is must be either numeric or character")
 }  else exprs<-x

if (!is.list(pathways) | is(pathways[[1]])!="pathway") stop("Argument 'pathways' is not a list of pathway-class objects")
if (convert) pathways<-convertIDs(pathways, IDs)
CheckNames(pathways, exprs)

if (type=="MA") out<-runTopologyGSAMulti(pathways, exprs, group, test=test, testCliques=testCliques, alpha=alpha, both.directions=both.directions, nperm=nperm, maxNodes = maxNodes)
 
if (type=="RNASeq") {
 if (is.null(norm.method)) {
  norm.method<-"TMM"
  cat("Normalization method was not specified. TMM used as default\n")
  out<-runTopologyGSA4RNASeq(pathways, exprs, group, test=test, testCliques=testCliques, alpha=alpha, both.directions=both.directions, nperm=nperm, norm.method, maxNodes = maxNodes)
  } else
  out<-runTopologyGSA4RNASeq(pathways, exprs, group, test=test, testCliques=testCliques, alpha=alpha, both.directions=both.directions, nperm=nperm, norm.method, maxNodes = maxNodes)
  }
return(out)
}

runTopologyGSAMulti<-function (pathways, exprs, group, test, testCliques, alpha, both.directions, nperm, maxNodes = maxNodes)
{
filterByTag<-function (tag, l) 
{
    isTagged <- sapply(l, function(x) x[[1]] == tag)
    lapply(l[isTagged], function(x) x[[2]])
}

    group<-factor(group)
    exp1<-t(exprs[,group==levels(group)[1]])
    exp2<-t(exprs[,group==levels(group)[2]])
    perm.num=nperm

    #if (!require(topologyGSA)) stop("library topologyGSA is missing")
    if (!any(test==c("mean","var"))) stop(paste("invalid test type:", test))
    #test <- switch(test, var = pathway.var.test, mean = pathway.mean.test)
    if (!is.null(maxNodes)) 
        pathways.small <- Filter(function(p) length(nodes(p)) <= maxNodes, pathways)
    f<-function(x){cat(x@title,"\n"); runTopologyGSADirected(x, test, exp1, exp2, alpha, both.directions, nperm)}
    log <- lapply(pathways.small, function(x) {
        tryCatch(list("ok", f(x)), error = function(e) list("err", 
            e))
    })
    out<-list(results = Filter(Negate(is.null), filterByTag("ok", 
        log)), errors = sapply(filterByTag("err", log), gettext))   
    
    out<-out[[1]]
    out<-data.frame(t(sapply(out, function(x) unlist(x[1:8]))))        
    if (testCliques) 
    {cat("Testing cliques..\n")
      cli<-lapply(pathways, function(path) {
      if (both.directions) path <- buildGraphNEL(nodes(path), edges(path), TRUE) else
          path <- buildGraphNEL(nodes(path), edges(path), FALSE)

      if (test=="mean") if (length(nodes(path))<=maxNodes) clique.mean.test(exp1, exp2, path, alpha, nperm)[c("p.value", "cliques")] else NULL
      if (test=="var")  if (length(nodes(path))<=maxNodes) clique.var.test(exp1, exp2, path, alpha)[c("p.value", "cliques")] else NULL
      })} else cli<-NULL
     out<-list(res=out, topo.sig=cli, degtable=NULL)
     class(out)<-c(class(out),  "topResultC","topResult")
    return(out)
}
#TBS
TBS<-function(x, group, pathways, type="MA", test=NULL, convert=TRUE, IDs="entrez",gene.stat="logFC", both.directions=TRUE, logFC.th=2, p.val.th=0.05, nperm=1000)
{
if (any(is(x) == "ExpressionSet")) {
exprs<-exprs(x) 
if (is.character(group) & ! any(group == colnames(pData(x)))) stop( paste( group,"is not present in phenoData"))  else
if (is.numeric(group)  & group > ncol(pData(x))) stop(paste("phenoData contain only", ncol(pData(x)), "columns")) else
if (is.numeric(group) | (is.character(group)) ) group<-pData(x)[,group] else 
 stop ("'group' is must be either numeric or character")
 }  else exprs<-x

if (!is.list(pathways) | is(pathways[[1]])!="pathway") stop("Argument 'pathways' is not a list of pathway-class objects")
if (convert) pathways<-convertIDs(pathways, IDs)
CheckNames(pathways, exprs)

if (type=="MA") out<-runTBS(pathways, exprs, group, gene.stat,both.directions,logFC.th=logFC.th, p.val.th=p.val.th, nperm=nperm)
if (type=="RNASeq") {
if (is.null(test)) {
  test<-"vstlimma"
  cat("test was not specified. 'vstlimma' used as default\n")
  out<-runTBS4RNASeq(pathways, exprs, group, gene.stat,both.directions,logFC.th=logFC.th, p.val.th=p.val.th, nperm=nperm, test=test)
  } else 
  out<-runTBS4RNASeq(pathways, exprs, group, gene.stat,both.directions,logFC.th=logFC.th, p.val.th=p.val.th, nperm=nperm, test=test)
  }
return(out)
}
#PWEA
PWEA<-function(x, group, pathways, type="MA", test=NULL, convert=TRUE, IDs="entrez", gene.stat="logFC", both.directions=TRUE, alpha=0.05, nperm=5000)
{
if (any(is(x) == "ExpressionSet")) {
exprs<-exprs(x) 
if (is.character(group) & ! any(group == colnames(pData(x)))) stop( paste( group,"is not present in phenoData"))  else
if (is.numeric(group)  & group > ncol(pData(x))) stop(paste("phenoData contain only", ncol(pData(x)), "columns")) else
if (is.numeric(group) | (is.character(group)) ) group<-pData(x)[,group] else 
 stop ("'group' is must be either numeric or character")
 }  else exprs<-x

if (!is.list(pathways) | is(pathways[[1]])!="pathway") stop("Argument 'pathways' is not a list of pathway-class objects")
if (convert) pathways<-convertIDs(pathways, IDs)
CheckNames(pathways, exprs)
if (type=="MA") out<-runPWEA(pathways, exprs, group, gene.stat, both.directions, alpha=alpha, nperm=nperm)
if (type=="RNASeq") {
if (is.null(test)) {
  test<-"vstlimma"
  cat("test was not specified. 'vstlimma' used as default\n")
  out<-runPWEA4RNASeq(pathways, exprs, group, gene.stat, both.directions, alpha=alpha, nperm=nperm, test=test)
  } else 
  out<-runPWEA4RNASeq(pathways, exprs, group, gene.stat, both.directions, alpha=alpha, nperm=nperm, test=test)
  }

return(out)
}
#Clipper
Clipper<-function(x, group, pathways, type="MA", convert=TRUE, IDs="entrez", both.directions=TRUE, test="mean", testCliques=FALSE, nperm=1000, norm.method=NULL)
{
if (any(is(x) == "ExpressionSet")) {
exprs<-exprs(x) 
if (is.character(group) & ! any(group == colnames(pData(x)))) stop( paste( group,"is not present in phenoData"))  else
if (is.numeric(group)  & group > ncol(pData(x))) stop(paste("phenoData contain only", ncol(pData(x)), "columns")) else
if (is.numeric(group) | (is.character(group)) ) group<-pData(x)[,group] else 
 stop ("'group' is must be either numeric or character")
 }  else exprs<-x

if (!is.list(pathways) | is(pathways[[1]])!="pathway") stop("Argument 'pathways' is not a list of pathway-class objects")
if (convert) pathways<-convertIDs(pathways, IDs)
CheckNames(pathways, exprs)
if (type=="MA") out<-runClipper(pathways, exprs, group, test, both.directions, testCliques=testCliques, nperm=nperm)
if (type=="RNASeq") {
 if (is.null(norm.method)) {
  norm.method<-"TMM"
  cat("Normalization method was not specified. TMM used as default\n")
  out<-runClipper4RNASeq(pathways, exprs, group,  method=test, both.directions, testCliques=testCliques, nperm=nperm, norm.method=norm.method)
  } else
  out<-runClipper4RNASeq(pathways, exprs, group,  method=test, both.directions, testCliques=testCliques, nperm=nperm, norm.method=norm.method) 
  }
return(out)
}