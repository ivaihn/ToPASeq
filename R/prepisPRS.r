#PRS

# Priprava permutacii
# all mena vsetkych
# de mena a fold-change DEG
# out matica nahodnych fold-change 
#preparePermsPRS<-function(all, de){
#ind<-as.numeric(all %in% names(de))
#perms.ind<-replicate(nperm, sample(ind))
#rownames(perms.ind)<-all
#perms<-apply(perms.ind, 2, function(x) {x[x==1]<-sample(de);x})
#return(perms)
#}

# PRS skore pre jednu drahu
PRSSingle<-function(path, de, all, perms){

 
weight<-PRSweights(path, de, all)
#g<-rownames(path)[rownames(path) %in% all]
#ind<-g %in% names(de)
#nf<-sum(ind)/length(g)
#expr<-ifelse(g %in% names(de), de[g], ifelse(g %in% all, 1, 0 ))

g<-rownames(path)
g<-lapply(g, function(x) strsplit(x, " ")[[1]])
g<-lapply(g, function(x) substr(x, regexpr(":",x)+1, nchar(x)))
ind<-sapply(g, function(x) any(as.character(x) %in% all))

indde<-sapply(g, function(x) any(as.character(x) %in% names(de)))
nf<-sum(indde)/length(g)
g<-g[ind ] 
expr<-sapply(g, function(x) if (any(as.character(x) %in% names(de))) max(de[as.character(x)], na.rm=TRUE) else if (any(as.character(x) %in% all)) 1 else 0)
obs<-sum(expr*weight)*nf

#random
gs<-unlist(g)[unlist(g) %in% all]
permsub<-perms[gs,, drop=FALSE]
gn<-unlist(sapply(seq_len(length(g)), function(i) rep(i, sum(g[[i]] %in% all))))
permsub<-apply(permsub, 2, function(x) tapply(x, gn, max, na.rm=T))


weight.rn<-apply(permsub,2, function(x) {

 weight<-setNames(rep(0, length(g)),rownames(path)[ind])
 if (length(g)>0 & length(g[x!=0])>=1)  weight[rownames(path)[ind][x!=0]]<-  downstreamCpp(path[ind,ind], rownames(path)[ind], rownames(path)[ind][x!=0])+1 #set, rownames(path)[ind1], rownames(path)[ind1][ind]
 weight[x==0]<-0
 return(weight) 
})

nf.rn<-colSums(permsub !=0)

rand<-colSums(weight.rn*permsub)*(nf.rn/length(g))

# normalization
obs<-(obs-mean(rand))/sd(rand)
rand<-(rand-mean(rand))/sd(rand)
p.value<-sum(rand >= obs)/length(rand)
res<-c(nPRS=obs, p.value=p.value)
return(res)
}
PRSweights<-function(path, de, all){

 g<-rownames(path)
 g<-lapply(g, function(x) strsplit(x, " ")[[1]])
 g<-lapply(g, function(x) substr(x, regexpr(":",x)+1, nchar(x)))

 ind1<-sapply(g, function(x) any(as.character(x) %in% all)) 
 g<-g[ind1]
 
 if (length(g)==0) stop("Pathway does not contain any measured genes")
 set<-path[ind1,ind1]
 
  ind<-sapply(g, function(x) any(as.character(x) %in% names(de)))  
# ind<-g %in% names(de)
 if (length(g)>0 & sum(ind)>=1) weight<-downstreamCpp(set, rownames(path)[ind1], rownames(path)[ind1][ind]) +1 else weight<-setNames(rep(0, length(g)),rownames(path)[ind1])
wei<-setNames(rep(0, length(g[ind1])),rownames(path)[ind1])
 wei[names(weight)]<-weight 
 return(wei)
 }
 
 
 
prs<-function(de, all, pathways, nperm){

perms<-preparePerms(all=all, de=de, nperm=nperm, method="PRS")

out<-catchErr(pathways, function(p) PRSSingle(p[[1]], de, all, perms))


if(length(out[[1]])>0){
out[[1]]<-data.frame(t(vapply(out[[1]], function(x) x, numeric(2))))
out[[1]]$q.value<-p.adjust(out[[1]]$p.value,"fdr") }

return(out)
}

collectWeightsPRS<-function(de, all, pathways){
out<-catchErr(pathways, function(p) PRSweights(p[[1]], de, all))
return(out[[1]])
}