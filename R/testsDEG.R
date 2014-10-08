testlimma<-function(exprs, group){
design <- model.matrix(~0+factor(group))
colnames(design) <- c("V1","V2")
fit <- lmFit(exprs, design)
contrast.matrix = makeContrasts(contrasts="V2-V1", levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)
fit2 = eBayes(fit2)
results <- decideTests(fit2)
deg.table = topTable(fit2, coef=1, adjust.method="BH", number=nrow(exprs), genelist=rownames(exprs))
return(deg.table)
}

testDESeq2<-function(exprs, group, betaPrior=TRUE){
#require(DESeq2)
ds <- DESeqDataSetFromMatrix(countData = exprs,
                             colData = data.frame(condition = factor(group)),
                             design = ~condition)
ds <- DESeq(ds, fitType = "parametric", test = "Wald",
           betaPrior = betaPrior)
DESeq2results <- results(ds, independentFiltering = FALSE,
                          cooksCutoff = FALSE)
DESeq2results <- DESeq2results[order(DESeq2results$padj), ]
deg.table<-as.data.frame(DESeq2results@listData)
rownames(deg.table)<-DESeq2results@rownames
deg.table<-data.frame(ID=rownames(deg.table), logFC=deg.table$log2FoldChange, t=deg.table$stat, pval=deg.table$pvalue, padj=deg.table$padj)
return(deg.table)
}


testDESeq2.perm<-function(exprs, group, nperm, betaPrior=TRUE){
#require(DESeq2)
ds <- DESeqDataSetFromMatrix(countData = exprs,
                             colData = data.frame(condition = factor(group)),
                             design = ~condition)
ds <- estimateSizeFactors(ds)
ds <- estimateDispersions(ds, fitType="parametric")
perm.test<-lapply(1:nperm, function(x) {
cat(x," ")
cD<-colData(ds)
cD[,1]<-sample(cD[,1])
colData(ds)<-cD
ds <- nbinomWaldTest(ds, betaPrior=betaPrior)
DESeq2results <- results(ds, independentFiltering = FALSE,
                          cooksCutoff = FALSE)
DESeq2results <- DESeq2results[order(DESeq2results$padj), ]
deg.table<-as.data.frame(DESeq2results@listData)
rownames(deg.table)<-DESeq2results@rownames
deg.table<-data.frame(ID=rownames(deg.table), logFC=deg.table$log2FoldChange, t=deg.table$stat, pval=deg.table$pvalue, padj=deg.table$padj)
return(deg.table)
})
return(perm.test)
}


#edgeR, limma
testvoomlimma<-function(exprs, group){
 nf <- edgeR::calcNormFactors.default(exprs, method = "TMM")
 design <- model.matrix(~group)
 v <- voom(exprs,design,plot=FALSE, lib.size = colSums(exprs) * nf)
 fit <- lmFit(v,design)
 fit <- eBayes(fit)
 deg.table<-topTable(fit,coef=2,number=nrow(exprs), adjust.method="BH",sort.by="p")
 deg.table<-data.frame(ID=rownames(deg.table), logFC=deg.table$logFC, t=deg.table$t, pval=deg.table$P.Value, padj=deg.table$adj.P.Val)
return(deg.table)
}

testvstlimma<-function(exprs, group){
 DESeq.cds = DESeq::newCountDataSet(countData = exprs, conditions = factor(group))
 DESeq.cds = DESeq::estimateSizeFactors(DESeq.cds)
 DESeq.cds = DESeq::estimateDispersions(DESeq.cds, method = "blind", fitType = "local")
 DESeq.vst = DESeq::getVarianceStabilizedData(DESeq.cds)
 DESeq.vst.fitlimma = lmFit(DESeq.vst, design = model.matrix(~factor(group)))
 DESeq.vst.fitbayes = eBayes(DESeq.vst.fitlimma)
 deg.table<-topTable(DESeq.vst.fitbayes,coef=2,number=nrow(exprs), adjust.method="BH",sort.by="p") 
 deg.table<-data.frame(ID=rownames(deg.table), 
                    logFC=deg.table$logFC, 
                        t=deg.table$t, 
                     pval=deg.table$P.Value,
                     padj=deg.table$adj.P.Val)
return(deg.table)
}

#edgeR
#dgelist <- DGEList(counts = count.matrix, group = classes)
#dgelist <- calcNormFactors(dgelist)
#dgelist <- estimateCommonDisp(dgelist)
#dgelist <- estimateTagwiseDisp(dgelist, trend = "movingave")
#exact.test = exactTest(dgelist)
#topTags(exact.test, n = nrow(count.matrix))