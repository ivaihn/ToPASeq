runPWEA4RNASeq <- function(sets,exprs, group, gene.stat="logFC", both.directions, alpha=0.05, nperm=5000, test="voomlimma"){
#require(genefilter)
cat("Preparing data..\n")    

if (test=="voomlimma") {
perm.test<-sapply(1:nperm, function(x) {
  cat(x," ")
  deg.table<-testvoomlimma(exprs, sample(group))
  test.all<-abs(deg.table$t)
  return(test.all)})
  
nf <- edgeR::calcNormFactors.default(exprs, method = 'TMM')
voom.data <- voom(exprs, design = model.matrix(~group),
                 lib.size = colSums(exprs) * nf,
                 plot = FALSE)
norm.exprs <- voom.data$E  

tif<-lapply(sets, function(x) {
if (suppressWarnings(CommonGenes(x, rownames(exprs)))) {
  return(NULL)
  }
  genes<-nodes(x)[nodes(x) %in% rownames(exprs)]
  TIF(norm.exprs[genes,],x, both.directions)})
tif.mean<-mean(unlist(tif))
tif.sd<-sd(unlist(tif))
 deg.table<-testvoomlimma(exprs, group)
obs.ttest<-abs(deg.table$t)
}


if (test=="vstlimma") {
perm.test<-sapply(1:nperm, function(x) {
  cat(x," ")
  deg.table<-testvstlimma(exprs, sample(group))
  test.all<-abs(deg.table$t)
  return(test.all)})
  
 DESeq.cds = DESeq::newCountDataSet(countData = exprs, conditions = factor(group))
 DESeq.cds = DESeq::estimateSizeFactors(DESeq.cds)
 DESeq.cds = DESeq::estimateDispersions(DESeq.cds, method = "blind", fitType = "local")
 DESeq.vst = DESeq::getVarianceStabilizedData(DESeq.cds)
 norm.exprs<-DESeq.vst
tif<-lapply(sets, function(x) {
if (suppressWarnings(CommonGenes(x, rownames(exprs)))) {
  return(NULL)
  }
  genes<-nodes(x)[nodes(x) %in% rownames(exprs)]
  TIF(norm.exprs[genes,],x, both.directions)})
tif.mean<-mean(unlist(tif))
tif.sd<-sd(unlist(tif))
 deg.table<-testvstlimma(exprs, group)
obs.ttest<-abs(deg.table$t)
}

if (test=="DESeq2") {
perm.test<-testDESeq2.perm(exprs, group, nperm)
perm.test<-sapply(perm.test, function(x) abs(x$t))

ds <- DESeqDataSetFromMatrix(countData = exprs,
                             colData = data.frame(condition = factor(group)),
                             design = ~condition)
ds <- DESeq(ds, fitType = "parametric", test = "Wald", betaPrior = TRUE)
vsd.matrix <- assay(varianceStabilizingTransformation(ds, blind=TRUE))
colnames(vsd.matrix) <- colnames(exprs)
norm.exprs<-vsd.matrix

tif<-lapply(sets, function(x) {
if (suppressWarnings(CommonGenes(x, rownames(exprs)))) {
  return(NULL)
  }
  genes<-nodes(x)[nodes(x) %in% rownames(exprs)]
  TIF(norm.exprs[genes,],x, both.directions)})
tif.mean<-mean(unlist(tif))
tif.sd<-sd(unlist(tif))
 deg.table<-testDESeq2(exprs, group)
obs.ttest<-abs(deg.table$t)
names(obs.ttest)<-rownames(exprs)
}
           
if (test=="limmaMA"){
cat("Preparing data..\n")
perm.test<-sapply(1:nperm, function(x) {
  if (x %% 100==0) cat(x," ") 
  if (x==nperm) cat("\n")
  test.all<-abs(testlimma(exprs, sample(group))$t)
  return(test.all)})

tif<-lapply(sets, function(x) {
if (suppressWarnings(CommonGenes(x, rownames(exprs)))) {
  return(NULL)
  }
  genes<-nodes(x)[nodes(x) %in% rownames(exprs)]
  TIF(exprs[genes,],x, both.directions)})
tif.mean<-mean(unlist(tif))
tif.sd<-sd(unlist(tif))
obs.ttest<-abs(testlimma(exprs, group)$t)
names(obs.ttest)<-rownames(exprs)
}                     
cat("Processing gene set:\n")                     
res<-sapply(sets, function(set) {
cat(set@title,"\n")
#observed
if (suppressWarnings(CommonGenes(set, rownames(exprs)))) {
  return(rep(NA,2))
  }
genes<-nodes(set)[nodes(set) %in% rownames(exprs)]
ind<-match(genes,rownames(exprs))
score.in<-obs.ttest[ind]^tif[[set@title]]
n.out<-nrow(exprs)-length(genes)#sum(!rownames(expr) %in% genes)
tif.out<-rnorm(n.out, tif.mean, tif.sd)
score.out<-obs.ttest[-ind]^tif.out
obs.es<-ES(score.in, score.out)
#permuted

perm.es<-apply(perm.test, 2, function(x) {
score.in<-x[ind]^tif[[set@title]]
score.out<-x[-ind]^tif.out
es<-ES(score.in, score.out)
return(es)
})
p<-mean(perm.es>=obs.es)

return(c(obs.es, p))
})
p.adj<-p.adjust(res[2,], "fdr")
res<-rbind(res,p.adj)
rownames(res)<-c("ES","p","p.adj")


tif<-sapply(tif, function(x) if (!is.null(x)) t(x))
if (gene.stat=="logFC") degt<-deg.table$logFC
if (gene.stat=="stats") degt<-deg.table$t
names(degt)<-rownames(exprs)
out<-list(res=data.frame(t(res)), topo.sig=tif, degtest=degt)
class(out)<-c(class(out), "topResultW","topResult")
return(out)
}
