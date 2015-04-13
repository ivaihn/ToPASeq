#############################
# Preparing Expression Data #
#############################

prepareData<-function(x, group=NULL, type, method, norm.method=NULL, test.method=NULL, nperm=NULL, ncores=NULL, p.th=0.05, logFC.th=1.5){

if (method %in% c("TAPPA", "DEGraph", "TopologyGSA", "clipper")) {
  if (type=="MA") { 
    out<-processMA(x, group) 
    } else
  if (type=="RNASeq") {
    if (is.null(norm.method)) norm.method<-"TMM"
    out<-processRNAseq(x, group, norm.method)
    } else stop("This inpout type is not supported for the selected method")  
  
  } else
if (method=="PWEA") {
  if (type=="RNAseq"){ 
   if (!is.null(norm.method)) warning("Normalization method ignored")
   if (is.null(test.method)) test.method<-"voomlimma"
   out<-testRNAseq(x, group, test.method)
   if (test.method=="DESeq2") perm<-testRNAseq(x, group, "DESeq2perm", nperm) else {
   test<-get(paste("test", test.method, sep=""))
   perm<-parDE(x, group,test , nperm, ncores)
   }
   out<-list(out, perm)
   }   else
  if (type=="MA"){
    out<-processMA(x, group)
    obs<-testMA(out[[1]], out[[2]])
    perm<-parDE(out[[1]], out[[2]], testMA, nperm, ncores)
    out<-list(obs, perm)
    } else
  if (type=="DEtable") {
  if (length(x)!=2) stop("The input data must be a list of length two (observed and random statistics)")
  if (!all(colnames(x[1])==c("ID", "logFC","t","pval","padj"))) stop('Colnames of the table differ from c("ID", "logFC","t","pval"."padj")"')
    out<-x
   } else stop("This inpout type is not supported for the selected method")  
  
  } else
if (method %in% c("PRS","SPIA")) {
  if (type=="RNASeq"){
   if (!is.null(norm.method)) warning("Normalization method ignored")
   if (is.null(test.method)) test.method<-"voomlimma"
   out<-testRNAseq(x, group, test.method)
   out<-applyTh(out, p.th, logFC.th)
  } else
  if (type=="MA"){
    out<-processMA(x, group)
    out<-testMA(out[[1]], out[[2]])
    out<-applyTh(out, p.th, logFC.th)
    if (length(out[[1]])==0) stop("Found NO differetially expressed genes") else message(paste("Found", length(out[[1]]), "differentially expressed genes"))
  } else
  if (type=="DEtable") {
    out<-applyTh(x, p.th, logFC.th)
  } else
  if (type=="DElist") {
    if (length(x)!=2) stop("The input data must be a list of length two (named statistics of differentially expressed genes and names of all genes)")
    checkDEandAll(x[[1]], x[[2]])
   } else stop("This inpout type is not supported for the selected method")  
  
} else stop(paste("Method",method,"is not implemented"))
return(out)
}

processMA<-function(x, group){
if (any(is(x) == "ExpressionSet")) {
exprs<-exprs(x) 
if (is.null(group)) stop("argument 'group' can not be NULL")
if (length(group)>1) stop("'group' must be of length one (number or name of the column of pData())")
if (is.character(group) & ! any(group == colnames(pData(x)))) stop( paste( group,"is not present in phenoData"))  else
if (is.numeric(group)  & group > ncol(pData(x))) stop(paste("phenoData contain only", ncol(pData(x)), "columns")) else
if (is.numeric(group) | (is.character(group)) ) group<-pData(x)[,group] else 
 stop ("'group' must be either numeric or character")
 } 

if (!any(is(x)%in% c("data.frame", "matrix"))) stop("Invalid 'x' type")

if (length(group)!=ncol(x)) stop("length of 'group' does not match ncol(x)")
if (nlevels(factor(group))!=2) stop(paste("'group' must contain two distict values, found", levels(factor(group))))
 
return(list(x, factor(group)))
}

processRNAseq<-function(x, group, norm.method){
if (! is.integer(x)) stop("Integer count data expected")
#normalizacia
if (norm.method=="TMM") {
nf = edgeR::calcNormFactors.default(x, method = 'TMM')
design <- model.matrix(~group)
voom.data = voom(x, design, lib.size = colSums(x) * nf, plot = FALSE)
norm.matrix = voom.data$E} else

if (norm.method=="DESeq2") {
 ds <- DESeqDataSetFromMatrix(countData = x,
                             colData = data.frame(condition = factor(group)),
                             design = ~condition)
ds <- DESeq(ds, fitType = "parametric", test = "Wald", betaPrior = TRUE)
vsd.matrix <- assay(varianceStabilizingTransformation(ds, blind=TRUE))
colnames(vsd.matrix) <- colnames(x)
norm.matrix<-vsd.matrix} else

if (norm.method=="rLog") {
 norm.matrix<-rlogTransformation(x)
} else 
if (norm.method=="none") norm.matrix<-x  else stop("Invalid normalization method")
return(list(x=norm.matrix, factor(group)))
}

testMA<-function(x, group){
cat(levels(group)[1],"denoted as 0 \n", levels(group)[2],"denoted as 1\n", "Contrasts: ", levels(group)[2], "-", levels(group)[1],"\n"    )
design <- model.matrix(~0+factor(group))
colnames(design) <- c("V1","V2")
fit <- limma::lmFit(x, design)
contrast.matrix = limma::makeContrasts(contrasts="V2-V1", levels=design)
fit2 = limma::contrasts.fit(fit, contrast.matrix)
fit2 = limma::eBayes(fit2)
results <- limma::decideTests(fit2)
deg.table = limma::topTable(fit2, coef=1, adjust.method="BH", number=nrow(x), genelist=rownames(x))
deg.table<-data.frame(ID=rownames(deg.table), logFC=deg.table$logFC, t=deg.table$t,  pval=deg.table$P.Value, padj=deg.table$adj.P.Val)
return(deg.table)
}

testRNAseq<-function(x, group, test.method, nperm=0){
if (test.method %in% c("DESeq2", "voomlimma", "vstlimma", "edgeR")) { 
deg.table<-get(paste("test", test.method, sep=""))(x,group)
} else 
if (test.method=="DESeq2perm") {
 deg.table<-testDESeq2perm(x, group, nperm)
 } else
stop("Invalid method for differential expression analysis")

return(deg.table)
}


parDE<-function(x, group, test, nperm, ncores=NULL){
if (is.null(ncores)) 
 ncores<-detectCores() else ncores=1
cl <- makeCluster(rep("localhost",ncores))
i<-1
ptest<-perm.test(i, test)
clusterExport(cl,deparse(substitute(test)))
tmp<-parSapply(cl, seq_len(nperm), get("ptest"), x, group)
stopCluster(cl)
return(tmp)
}

perm.test<-function(i, test){
function(i,exprs, group){test(exprs, sample(group))}
}


applyTh<-function(x, p.th, logFC.th){

if (!is.na(p.th) & !is.na(logFC.th)) select<-(x$pval < p.th) & (abs(x$logFC) > logFC.th) else 
if (!is.na(p.th)) select<- x$pval < p.th else
if (!is.na(logFC.th)) select<-abs(x$logFC) > logFC.th else select<-rep(FALSE, nrow(x))

de<-setNames(x$logFC[select], x$ID[select])
out<-list(de=de, all=as.character(x$ID))
return(out)
}



 

######################
# Preparing pathways #
######################
preparePathways<-function(pathways, method, both.directions, genes, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="entrez", convertBy=NULL, EdgeAttrs=NULL){

#konverzia indetifikatorov
if (is.null(convertBy)) {
 if (any(convertTo==c("entrez","symbol"))) pathways<-sapply(pathways, convertIdentifiers, convertTo) else
 if (convertTo!="none") stop("Invalid conversion specified") 
 } else pathways<-sapply(pathways, convertIdentifiersByVector, convertBy, "user_defined") 

CheckNames(pathways, genes)

N<-length(pathways)
# privelke drahy
if (!is.null(maxNodes)) pathways<-BigPaths(pathways, maxNodes) 
if (!is.null(commonTh)) pathways<-CommonGenes(pathways, genes, commonTh)
if (!is.null(minEdges)) pathways<-FewEdges(pathways, minEdges)


  pathways<-lapply(pathways, function(path) list(transformPathway(path, method=method, both.directions, EdgeAttrs), path@title))

if (filterSPIA) pathways<-filterSPIA(pathways)
message(paste(N-length(pathways), " pathways were filtered out"))

return(pathways)
}

# pomocne
preparePerms<-function(de=NULL, all=NULL, x=NULL, group=NULL, nperm=NULL, test=NULL, method){
if (method=="PRS"){
  print(head(all))
  print(head(de))
    ind<-as.numeric(all %in% names(de))
    perms.ind<-replicate(nperm, sample(ind))
    rownames(perms.ind)<-all
    perms<-apply(perms.ind, 2, function(x) {x[x==1]<-sample(de);x})
 }
if (method=="PWEA"){
  perms<-replicate(nperm, test(x, sample(group))$stat)
 }
if (method=="TopologyGSA") {
  perms<-replicate(nperm, sample(group))
 } 
 return(perms)
}

CheckNames<-function(pathway, expr){

IDmatchsum<-sapply(pathway, function(x) sum(nodes(x) %in% expr))
IDmatchmean<-sapply(pathway, function(x) mean(nodes(x) %in% expr))
if (sum(IDmatchsum)==0) stop("Gene labels and node labels do not match. Please, correct your gene identifiers\n")
cat(sum(IDmatchsum),"node labels mapped to the expression data\n")
cat("Average coverage", mean(IDmatchmean,na.rm=TRUE)*100,"%\n")
cat(sum(IDmatchsum==0)," (out of ",length(pathway),") pathways without a mapped node\n", sep="")
}

transformPathway<-function(x, method, both.directions=TRUE, EdgeAttrs=NULL){
if (! any(class(x)=="pathway")) stop("x must be an object of 'pathway'-class")
if (is.null(EdgeAttrs)) EdgeAttrs<-makeDefaultEdgeData() 


if (method=="TAPPA") {
  x<-as(x,"graphNEL") 
  x<-as(x,"matrix")
  x<-x+t(x)
  x[x>1]<-1
  x[lower.tri(x)]<-0
  diag(x)<-1
  }
if (method=="PRS") {
  if (both.directions) x<-buildGraphNEL(nodes(x), edges(x), TRUE) else
  x<-buildGraphNEL(nodes(x), edges(x), FALSE)
  x<-as(x,"matrix")
}
if (any(method==c("PWEA","TopologyGSA","clipper","DEGraphNoSigns"))) {
  if (both.directions) x<-buildGraphNEL(nodes(x), edges(x), TRUE) else
  x<-buildGraphNEL(nodes(x), edges(x), FALSE)
}

if (method=="DEGraph") {
  if (both.directions) g<-buildGraphNEL(nodes(x), edges(x), TRUE) else
  g<-buildGraphNEL(nodes(x), edges(x), FALSE)
  eA<-merge(EdgeAttrs[[1]], EdgeAttrs[[2]], by.x=2, by.y=1, all=TRUE)
  pos<-as.character(eA[,2][!is.na(eA[,2]) & eA[,"beta"]==1])
  neg<-as.character(eA[,2][!is.na(eA[,2]) & eA[,"beta"]==-1])
  neu<-as.character(eA[,2][!is.na(eA[,2]) & eA[,"beta"]==0])

  signMat<-as(g,"matrix")
  signMat[,]<-0

  e<-edges(x)
  posedg<-e[,1:2][e[,4] %in% pos, ]
  posind<-nrow(signMat)*(match(posedg[,2], colnames(signMat))-1)+match(posedg[,1], colnames(signMat))
  signMat[posind] <-  1
  negedg<-e[,1:2][e[,4] %in% neg, ]
  negind<-nrow(signMat)*(match(negedg[,2], colnames(signMat))-1)+match(negedg[,1], colnames(signMat))
  signMat[negind] <- -1

  neuedg<-e[,1:2][e[,4] %in% neu, ]
  g<-removeEdge(from=neuedg[,1], to=neuedg[,2], g)
  g@graphData$signMat<-signMat
  x<-g
}

if (method=="SPIA") { 
 x<-prepareSPIA2(x, both.directions, EdgeAttrs)
 x<-getdatp(x, EdgeAttrs$beta$rel, EdgeAttrs$beta$beta)
}
return(x)

}

# priprava topologii pre SPIA, POZOR bez konverzie identifikatorov
prepareSPIA2<-function (p, both.directions, edgeAttrs=makeDefaultEdgeData()){
    
        es <- graphite::edges(p)
        ns <- graphite::nodes(p)
        spiaEdges<-edgeAttrs[[1]]
        if (!all(levels(es[,4]) %in% spiaEdges[,1] )) stop("Unexpected edge type ", levels(es[,4])[!levels(es[,4]) %in% spiaEdges[,1]], " Please see '?makeDefaultEdgeData' for explanation\n")       
        es <- merge(es, spiaEdges, all.x = TRUE)[c("src", "dest", "direction", "spiaType")]
        l <- sapply(edgeAttrs[[2]][,1], simplify = FALSE, USE.NAMES = TRUE,
            function(edgeType) {
                est <- es[es[, 4] == edgeType, , drop = FALSE]
                if (both.directions) gnl <- buildGraphNEL(ns, est, TRUE) else
                 gnl<-buildGraphNEL(ns, est, FALSE)
                t(as(gnl, "matrix"))
            })
        l$title <- p@title
        l$nodes <- ns
        l$NumberOfReactions <- 0
        
        return(l)
    
}

filterSPIA<-function(pathways){
    pathways <- Filter(function(p) sum(abs(p))==0, pathways)

    return(pathways)    
    
   }
   
FewEdges<-function (pathways, minEdges) {
    pathways <- Filter(function(p) nrow(edges(p)) > minEdges , pathways)
    return(pathways)
} 
   
BigPaths<-function (pathways, maxNodes) {
    pathways <- Filter(function(p) length(nodes(p)) <= maxNodes, pathways)
    return(pathways)
}

CommonGenes<-function (pathways, genes, threshold) {
    pathways<-Filter(function(p) length(intersect(nodes(p), genes)) >= threshold, pathways)
    return(pathways)
}

catchErr<-function (l, f){
    log <- lapply(l, function(x) {tryCatch(list("ok", f(x)), error = function(e) list("err", e))})
    list(results = Filter(Negate(is.null), filterByTag("ok", log)), errors = sapply(filterByTag("err", log), gettext))
}

filterByTag<-function (tag, l) 
{
    isTagged <- sapply(l, function(x) x[[1]] == tag)
    lapply(l[isTagged], function(x) x[[2]])
}



    checkPkgVersion<-function (name, min_version)
{
    version <- package_version(installed.packages()[name, "Version"])
    if (version < package_version(min_version))
        stop("the installed ", name, " version is too old (need at least ",
            min_version, ")")
}

   
###########
# TAPPA
#deg.table<-testlimma(exprs, group)
#if (gene.stat=="logFC") degt<-deg.table$logFC
#if (gene.stat=="stats") degt<-deg.table$t
#names(degt)<-deg.table$ID
#out<-list(res=t(PCI), topo.sig=NULL, degtest=degt)
#class(out)<-c(class(out), "topResultE","topResult")