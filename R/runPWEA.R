runPWEA <-
function(sets, expr, group, gene.stat="logFC", both.directions, alpha=0.05, nperm=5000){

group<-factor(group)
if (!nlevels(group)==2) stop("Group vector has not two levels")

cat(levels(group)[1],"denoted as 0 \n", levels(group)[2],"denoted as 1\n", "Contrasts: ", levels(group)[1], "-", levels(group)[2],"\n"    )

#testlimma(expr, group)
#genefilter
cat("Preparing data..\n")
perm.test<-sapply(1:nperm, function(x) {
  if (x %% 100==0) cat(x," ") 
  if (x==nperm) cat("\n")
  test.all<-abs(testlimma(expr, sample(group))$t)
  return(test.all)})

tif<-lapply(sets, function(x) {
if (suppressWarnings(CommonGenes(x, rownames(expr)))) {
  return(NULL)
  }
  genes<-nodes(x)[nodes(x) %in% rownames(expr)]
  TIF(expr[genes,],x, both.directions)})
tif.mean<-mean(unlist(tif))
tif.sd<-sd(unlist(tif))
obs.ttest<-abs(testlimma(expr, group)$t)
names(obs.ttest)<-rownames(expr)
                     
cat("Processing gene set:\n")                     
res<-sapply(sets, function(set) {
cat(set@title,"\n")
#observed
if (suppressWarnings(CommonGenes(set, rownames(expr)))) {
  return(rep(NA,2))
  }
genes<-nodes(set)[nodes(set) %in% rownames(expr)]
#obs.test<-abs(rowttests(expr, as.factor(gr))$statistic)
ind<-match(genes,rownames(expr))
score.in<-obs.ttest[ind]^tif[[set@title]]
n.out<-nrow(expr)-length(genes)#sum(!rownames(expr) %in% genes)
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

tif<-sapply(tif, function(x) if (is.matrix(x)) t(x))
deg.table<-testlimma(expr, group)
if (gene.stat=="logFC") degt<-deg.table$logFC
if (gene.stat=="stats") degt<-deg.table$t
names(degt)<-rownames(expr)
out<-list(res=data.frame(t(res)), topo.sig=tif, degtest=degt)
class(out)<-c(class(out), "topResultW","topResult")
return(out)
}
