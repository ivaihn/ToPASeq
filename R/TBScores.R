#TBScores <-
#function(set, nodes.info, de, all, nperm=1000){
#TBS.obs<-sum(nodes.info[1,]*nodes.info[2,], na.rm=TRUE)
#de.path<-length(de[names(de) %in% nodes(set) ])
#TBS.obs<-TBS.obs * de.path/length(all)
#
#TBS.perm<-vapply(1:nperm, function(x) {
#de.rn<-de
#names(de.rn)<-sample(all, length(de))
#
#ni.rn<-processNodes2(set, de.rn, all, down.list)
#de.path.rn<-length(de.rn[names(de.rn) %in% nodes(set) ])
#TBS.rn<-sum(ni.rn[1,]*ni.rn[2,]) * length(de.path.rn)/length(all)
#return(TBS.rn)
#}, 0)
#
#
#nTBS.obs<-(TBS.obs-mean(TBS.perm))/sd(TBS.perm)
#nTBS.perm<-(TBS.perm-mean(TBS.perm))/sd(TBS.perm)
#pval<-mean(nTBS.obs<nTBS.perm)
#
#return(list(normalized.observed=nTBS.obs, normalized.random=nTBS.perm, observed=TBS.obs, random=TBS.perm, p.value=pval, nodes.info=nodes.info))
#}

TBScores<-function(sets, de, all, down.list){
tbs<-sapply(sets, function(set) {
 nodes.info<-processNodes(set, de, all, down.list)
 sc<-sum(nodes.info[1,]*nodes.info[2,], na.rm=TRUE)
 nde<-length(de[names(de) %in% nodes(set) ])
 res<-sc * nde/length(all) 
 return(res)
 })
return(tbs) 
}
