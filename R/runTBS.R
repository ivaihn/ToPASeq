runTBS <-
function(sets, exprs, group,gene.stat="logFC", both.directions, logFC.th=2, p.val.th=0.05, nperm=1000 ) {
#if (!require(limma)) stop("Please, install limma package")
group<-factor(group)
if (!nlevels(group)==2) stop("Group vector has not two levels")

cat(levels(group)[1],"denoted as 0 \n", levels(group)[2],"denoted as 1\n", "Contrasts: ", levels(group)[1], "-", levels(group)[2],"\n"    )

deg.table<-testlimma(exprs, group)

de<-deg.table[abs(deg.table$logFC)>logFC.th & deg.table$P.Value<p.val.th, "logFC"]
names(de)<-deg.table[abs(deg.table$logFC)>logFC.th & deg.table$P.Value<p.val.th, "ID"]
cat("Found", length(de), "differentially expressed genes\n")
all<-deg.table$ID
cat("Preparing permutation table and downstream list\n")
perm.table<-sapply(1:nperm, function(x) sample(all, length(de)))
down.list<-sapply(sets, function(set) sapply(nodes(set), function(x) downstream.nodes(x,set,nodes(set), both.directions)))
names(down.list)<-sapply(sets, function(set) set@title)

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
p<-sapply(1:length(TBS.obs), function(i) sum(TBS.obs.norm[i]<=TBS.rn.norm[i,])/nperm)
p.adj<-p.adjust(p, "fdr")
TBS<-data.frame(TBS.obs.norm, p, p.adj)

ni<-sapply(sets, function(set) {
 nodes.info<-processNodes(set, de, all, down.list)[c(2,1),]
 } )
if (gene.stat=="logFC") degt<-deg.table$logFC
if (gene.stat=="stats") degt<-deg.table$t
names(degt)<-deg.table$ID
degt[! (abs(deg.table$logFC)>logFC.th & deg.table$P.Value<p.val.th)]<-1
out<-list(res=data.frame(TBS), topo.sig=ni, degtest=degt)
class(out)<-c(class(out), "topResultW","topResult")
return(out)
}
