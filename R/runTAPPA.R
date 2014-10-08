runTAPPA <-
function(sets, exprs, group, gene.stat="logFC", both.directions, normalize=TRUE, verbose=FALSE){   
group<-factor(group)
cat(levels(group)[1], "denoted as 1 ", "\n", levels(group)[2], "denoted as 2 \n")
group<-as.numeric(group)
if (normalize) {
  norm.data<-normalizeTAPPA(exprs)
  } else {
  norm.data<-exprs
  }

PCI<-sapply(sets, function(x) {
if (verbose) cat(x@title,"\n")
if (suppressWarnings(CommonGenes(x, rownames(exprs)))) {
  return(rep(NA,9))
  }
nodes<-x@nodes[x@nodes %in% rownames(norm.data)]
scores1<-PCIscores(norm.data[nodes,group==1], x, both.directions)
scores2<-PCIscores(norm.data[nodes,group==2], x, both.directions)
if (length(scores1)>2 & length(scores2)>2) return(c(sum(!is.na(scores1)), median(scores1), range(scores1), sum(!is.na(scores2)), median(scores2), range(scores2),wilcox.test(scores1,scores2)$p.value)) else return(c(rep(NA,9)))
})

rownames(PCI)<-c("valid1","median1", "min1", "max1","valid2", "median2", "min2", "max2", "p")
p.adj<-p.adjust(PCI["p",],"fdr")
PCI<-rbind(PCI,p.adj)

deg.table<-testlimma(exprs, group)
if (gene.stat=="logFC") degt<-deg.table$logFC
if (gene.stat=="stats") degt<-deg.table$t
names(degt)<-deg.table$ID
out<-list(res=t(PCI), topo.sig=NULL, degtest=degt)
class(out)<-c(class(out), "topResultE","topResult")
return(out)
}


