runTBS4RNASeq <- function(sets, exprs, group, gene.stat="logFC",both.directions, logFC.th=2, p.val.th=0.05, nperm=1000, test=c("vstlimma","DESeq2")) {

if (test=="DESeq2") deg.table<-testDESeq2(exprs, group)
if (test=="voomlimma") deg.table<-testvoomlimma(exprs, group)
if (test=="vstlimma") deg.table<-testvstlimma(exprs, group)
if (test=="limmaMA") deg.table<-testlimma(exprs,group)
if (!test %in% c("voomlimma","vstlimma", "DESeq2", "limmaMA")) stop('Invalid "test" argument. Use either "voomlimma","vstlimma" ,"limmaMA" or "DESeq2"')

de<-deg.table[abs(deg.table[,2])>logFC.th & deg.table[,5]<p.val.th, 2]
names(de)<-rownames(deg.table)[abs(deg.table[,2])>logFC.th & deg.table[,5]<p.val.th]
cat("Found", length(de), "differentially expressed genes\n")
all<-rownames(deg.table)
     

cat("Preparing permutation table and downstream list\n")
perm.table<-sapply(1:nperm, function(x) sample(all, length(de)))
down.list<-sapply(sets, function(set) sapply(nodes(set), function(x) downstream.nodes(x,set,nodes(set), both.directions)))

cat("Observed scores..\n")
TBS.obs<-TBScores(sets, de, all, down.list)
cat("Random scores..\n")
TBS.rn<-sapply(1:nperm, function(i) {
if (i %% 100==0) cat(i," ") 
if (i==nperm) cat("\n")
de.rn<-setNames(de, perm.table[,i])
TBScores(sets, de.rn, all, down.list)
})

cat("Normalization and p-values...\n")
TBS.obs.norm<-sapply(1:length(TBS.obs), function(i) (TBS.obs[i]-mean(TBS.rn[i,]))/sd(TBS.rn[i,]))
TBS.rn.norm<-t(sapply(1:length(TBS.obs), function(i) (TBS.rn[i,]-mean(TBS.rn[i,]))/sd(TBS.rn[i,])))

p<-sapply(1:length(TBS.obs), function(i) sum(TBS.obs.norm[i]>=TBS.rn.norm[i,])/nperm)
p.adj<-p.adjust(p, "fdr")
TBS<-data.frame(TBS.obs.norm, p, p.adj)

ni<-sapply(sets, function(set) {
 nodes.info<-processNodes(set, de, all, down.list)[c(2,1),]
 }             )
if (gene.stat=="logFC") degt<-deg.table$logFC
if (gene.stat=="stats") degt<-deg.table$t
names(degt)<-deg.table$ID
degt[! (abs(deg.table$logFC)>logFC.th & deg.table$P.Value<p.val.th)]<-1
out<-list(res=data.frame(TBS), topo.sig=ni, degtest=degt)
class(out)<-c(class(out), "topResultW","topResult")
return(out)
}
