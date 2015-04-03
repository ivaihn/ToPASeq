testDESeq2<-function(x, group){
ds <- DESeq2::DESeqDataSetFromMatrix(countData = x,
                             colData = data.frame(condition = factor(group)),
                             design = ~condition)
 ds <- DESeq2::DESeq(ds, fitType = "parametric", test = "Wald", betaPrior = TRUE)
 DESeq2results <- DESeq2::results(ds, independentFiltering = FALSE, cooksCutoff = FALSE)
 DESeq2results <- DESeq2results[order(DESeq2results$padj), ]
 deg.table<-as.data.frame(DESeq2results@listData)
 rownames(deg.table)<-DESeq2results@rownames
 deg.table<-data.frame(ID=rownames(deg.table), logFC=deg.table$log2FoldChange, t=deg.table$stat, pval=deg.table$pvalue, padj=deg.table$padj)
 return(deg.table)
 }
 
testDESeq2perm<-function(x, group, nperm){
 ds <- DESeq2::DESeqDataSetFromMatrix(countData = x, colData = data.frame(condition = factor(group)), design = ~condition)
 ds <- DESeq2::estimateSizeFactors(ds)
 ds <- DESeq2::estimateDispersions(ds, fitType="parametric")
 perm.test<-replicate(nperm, {
 cD<-colData(ds)
 cD[,1]<-sample(cD[,1])
 colData(ds)<-cD
 ds <- DESeq2::nbinomWaldTest(ds, betaPrior=TRUE)
 DESeq2results <- DESeq2::results(ds, independentFiltering = FALSE, cooksCutoff = FALSE)
 DESeq2results <-DESeq2results[order(DESeq2results$padj), ]
 deg.table<-as.data.frame(DESeq2results@listData)
 rownames(deg.table)<-DESeq2results@rownames
 deg.table<-data.frame(ID=rownames(deg.table), logFC=deg.table$log2FoldChange, t=deg.table$stat, pval=deg.table$pvalue, padj=deg.table$padj)
 return(deg.table)
 })
 return(perm.test)
 }
 
testvoomlimma<-function(x, group){
nf = edgeR::calcNormFactors.default(x, method = 'TMM')
 design <- model.matrix(~group)
 voom.data = limma::voom(x, design, lib.size = colSums(x) * nf, plot = FALSE)
 fit <- limma::lmFit(x,design)
 fit <- limma::eBayes(fit)
 deg.table<-limma::topTable(fit,coef=2,number=nrow(x), adjust.method="BH",sort.by="p")
 deg.table<-data.frame(ID=rownames(deg.table), logFC=deg.table$logFC, t=deg.table$t, pval=deg.table$P.Value, padj=deg.table$adj.P.Val)
 return(deg.table)
 }

testvstlimma<-function(x, group){
DESeq.cds = DESeq::newCountDataSet(countData = x, conditions = factor(group))
 DESeq.cds = DESeq::estimateSizeFactors(DESeq.cds)
 DESeq.cds = DESeq::estimateDispersions(DESeq.cds, method = "blind", fitType = "local")
 DESeq.vst = DESeq::getVarianceStabilizedData(DESeq.cds)
 DESeq.vst.fitlimma = limma::lmFit(DESeq.vst, design = model.matrix(~factor(group)))
 DESeq.vst.fitbayes = limma::eBayes(DESeq.vst.fitlimma)
 deg.table<-limma::topTable(DESeq.vst.fitbayes,coef=2,number=nrow(x), adjust.method="BH",sort.by="p") 
 deg.table<-data.frame(ID=rownames(deg.table), logFC=deg.table$logFC, t=deg.table$t, pval=deg.table$P.Value, padj=deg.table$adj.P.Val)
return(deg.table) 
}

testedgeR<-function(x, group){
 dgelist <- edgeR::DGEList(counts = x, group = group)
 dgelist <- edgeR::calcNormFactors(dgelist)
 dgelist <- edgeR::estimateCommonDisp(dgelist)
 dgelist <- edgeR::estimateTagwiseDisp(dgelist, trend = "movingave")
 exact.test <- edgeR::exactTest(dgelist)
 deg.table<-edgeR::topTags(exact.test, n = nrow(x))[[1]]
 deg.table<-data.frame(ID=rownames(deg.table), logFC=deg.table$logFC, t=NA , pval=deg.table$PValue, padj=deg.table$FDR)
 return(deg.table)
 }