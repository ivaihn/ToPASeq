runTAPPA4RNASeq<-function(sets, count.matrix, group, gene.stat="logFC",both.directions, normalize=TRUE, verbose=FALSE, norm.method=c("TMM","DESeq2", "none")){
count.matrix <- count.matrix[apply(count.matrix, 1, sum) > 0, ]

if (norm.method[1]=="TMM") {
nf = edgeR::calcNormFactors.default(count.matrix, method = 'TMM')
voom.data = voom(count.matrix, design = model.matrix(~group),
                 lib.size = colSums(count.matrix) * nf,
                 plot = FALSE)
norm.matrix = voom.data$E}

if (norm.method[1]=="DESeq") {
 ds <- DESeqDataSetFromMatrix(countData = count.matrix,
                             colData = data.frame(condition = factor(group)),
                             design = ~condition)
ds <- DESeq(ds, fitType = "parametric", test = "Wald", betaPrior = TRUE)
vsd.matrix <- assay(varianceStabilizingTransformation(ds, blind=TRUE))
colnames(vsd.matrix) <- colnames(count.matrix)
norm.matrix<-vsd.matrix}

if (norm.method[1]=="none") norm.matrix<-count.matrix
out<-runTAPPA(sets, norm.matrix, group, gene.stat, both.directions, normalize, verbose)
return(out)
}


runDEGraph4RNASeq<-function(sets, count.matrix, group, gene.stat="logFC", overall="biggest",  norm.method=c("TMM","DESeq2", "none"), useInteractionSigns, both.directions){
count.matrix <- count.matrix[apply(count.matrix, 1, sum) > 0, ]

if (norm.method[1]=="TMM") {
nf = edgeR::calcNormFactors.default(count.matrix, method = 'TMM')
voom.data = voom(count.matrix, design = model.matrix(~group),
                 lib.size = colSums(count.matrix) * nf,
                 plot = FALSE)
norm.matrix = voom.data$E}

if (norm.method[1]=="DESeq2") {
 ds <- DESeqDataSetFromMatrix(countData = count.matrix,
                             colData = data.frame(condition = factor(group)),
                             design = ~condition)
ds <- DESeq(ds, fitType = "parametric", test = "Wald", betaPrior = TRUE)
vsd.matrix <- assay(varianceStabilizingTransformation(ds, blind=TRUE))
colnames(vsd.matrix) <- colnames(count.matrix)
norm.matrix<-vsd.matrix}

if (norm.method[1]=="none") norm.matrix<-count.matrix

group<-factor(group)
res<-runDEGraphSigned(sets, norm.matrix, group, useInteractionSigns, both.directions, maxNodes = 150)
res<-processDEGraph(res[[1]], overall)

if (norm.method[1]=="TMM") deg.table<-testvoomlimma(count.matrix, group)
if (norm.method[1]=="DESeq2") deg.table<-testDESeq2(count.matrix, group)
if (norm.method[1]=="none") deg.table<-testDESeq2(count.matrix, group)


  if (gene.stat=="logFC") degt<-deg.table$logFC
  if (gene.stat=="stats") degt<-deg.table$t
  names(degt)<-deg.table$ID
  out<-list(res=res, topo.sig=NULL, degtest=degt)
  return(out)
}


runTopologyGSA4RNASeq<-function(sets, count.matrix, group, test="mean", testCliques=FALSE, alpha=0.05, both.directions, nperm=10000, norm.method=c("TMM","DESeq2", "none"), maxNodes=maxNodes){


if (norm.method[1]=="TMM") {
nf = edgeR::calcNormFactors.default(count.matrix, method = 'TMM')
voom.data = voom(count.matrix, design = model.matrix(~group),
                 lib.size = colSums(count.matrix) * nf,
                 plot = FALSE)
norm.matrix = voom.data$E}

if (norm.method[1]=="DESeq") {
 ds <- DESeqDataSetFromMatrix(countData = count.matrix,
                             colData = data.frame(condition = factor(group)),
                             design = ~condition)
ds <- DESeq(ds, fitType = "parametric", test = "Wald", betaPrior = TRUE)
vsd.matrix <- assay(varianceStabilizingTransformation(ds, blind=TRUE))
colnames(vsd.matrix) <- colnames(count.matrix)
norm.matrix<-vsd.matrix}

if (norm.method[1]=="none") norm.matrix<-count.matrix

group<-factor(group)
out<-runTopologyGSAMulti(sets, norm.matrix, group, test=test, testCliques=testCliques, alpha=alpha, both.directions=both.directions, nperm=nperm, maxNodes=maxNodes )
return(out)
}


runClipper4RNASeq<-function(sets, count.matrix, group,  method="mean", both.directions, testCliques=FALSE, nperm=NULL, norm.method=c("TMM","DESeq2", "none")){
count.matrix <- count.matrix[apply(count.matrix, 1, sum) > 0, ]

if (norm.method[1]=="TMM") {
nf = edgeR::calcNormFactors.default(count.matrix, method = 'TMM')
voom.data = voom(count.matrix, design = model.matrix(~group),
                 lib.size = colSums(count.matrix) * nf,
                 plot = FALSE)
norm.matrix = voom.data$E}

if (norm.method[1]=="DESeq2") {
 ds <- DESeqDataSetFromMatrix(countData = count.matrix,
                             colData = data.frame(condition = factor(group)),
                             design = ~condition)
ds <- DESeq(ds, fitType = "parametric", test = "Wald", betaPrior = TRUE)
vsd.matrix <- assay(varianceStabilizingTransformation(ds, blind=TRUE))
colnames(vsd.matrix) <- colnames(count.matrix)
norm.matrix<-vsd.matrix}

if (norm.method[1]=="none") norm.matrix<-count.matrix

group<-factor(as.numeric(factor(group)))
out<-runClipper(sets, norm.matrix, group, method, both.directions, testCliques, nperm )
return(out)
}


runSPIA4RNASeq<-function (pathways, exprs, group, gene.stat="logFC", both.directions, convert=TRUE, logFC.th=2, p.val.th=0.05, IDs, test, edgeAttrs, ...)
{

path.info<-prepareSPIA2(pathways,"pathways", both.directions=both.directions, convert=convert, IDs=IDs, edgeAttrs=edgeAttrs)

#if (!require(limma)) stop("Please, install limma package")
group<-factor(group)
if (!nlevels(group)==2) stop("Group vector has not two levels")
cat(levels(group)[1],"denoted as 0 \n", levels(group)[2],"denoted as 1\n", "Contrasts: ", levels(group)[2], "-", levels(group)[1],"\n"    )

if (test=="DESeq2") deg.table<-testDESeq2(exprs, group)
if (test=="voomlimma") deg.table<-testvoomlimma(exprs, group)
if (test=="vstlimma") deg.table<-testvstlimma(exprs, group)
if (test=="limmaMA") deg.table<-testlimma(exprs,group)
if (!test %in% c("voomlimma","vstlimma", "DESeq2", "limmaMA")) stop('Invalid "test" argument. Use either "voomlimma","vstlimma" ,"limmaMA" or "DESeq2"')


de<-       deg.table[abs(deg.table$logFC)>logFC.th & deg.table$pval<p.val.th, "logFC"]
names(de)<-deg.table[abs(deg.table$logFC)>logFC.th & deg.table$pval<p.val.th, "ID"]
cat("Found", length(de), "differentially expressed genes\n")

all<-deg.table$ID
    res<-graphite.SPIA(de, all, "pathways", beta=setNames(edgeAttrs[[2]]$beta, edgeAttrs[[2]][,1]))
     names(res)[7:9]<-c("p","pFdr","pFWER")
    rownames(res)<-res[,1]
    res<-res[,-1]
    if (gene.stat=="logFC") degt<-deg.table$logFC
    if (gene.stat=="stats") degt<-deg.table$t
    names(degt)<-deg.table$ID
    
    #path.info<-prepareSPIA2(pathways,"pathways", convert=convert, IDs=IDs, reduce=FALSE, edgeAttrs=edgeAttrs) 
    tsig<-spiaPF(de,all,"pathways")
    out<-list(res=res, topo.sig=tsig, degtest=degt)
    class(out)<-c(class(out), "topResultW","topResult")
return(out)
}