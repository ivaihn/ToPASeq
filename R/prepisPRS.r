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
g<-rownames(path)[rownames(path) %in% all]
ind<-g %in% names(de)
nf<-sum(ind)/length(g)

expr<-ifelse(g %in% names(de), de[g], ifelse(g %in% all, 1, 0 ))
obs<-sum(expr*weight)*nf

#random

weight.rn<-apply(perms[g,, drop=FALSE],2, function(x) {

 weight<-setNames(rep(0, length(g)),g)
 if (length(g)>0 & length(g[x!=0])>=1)  weight[g[x!=0]]<-  downstreamCpp(path[g,g], g, g[x!=0])+1
 weight[x==0]<-0
 return(weight) 
})

nf.rn<-colSums(perms[g,,drop=FALSE] !=0)

rand<-colSums(weight.rn*perms[g,])*(nf.rn/length(g))

# normalization
obs<-(obs-mean(rand))/sd(rand)
rand<-(rand-mean(rand))/sd(rand)
p.value<-sum(rand >= obs)/length(rand)
res<-c(nPRS=obs, p.value=p.value)
return(res)
}

PRSweights<-function(path, de, all){
 
 g<-rownames(path)[rownames(path) %in% all]
 if (length(g)==0) stop("Pathway does not contain any measured genes")
 set<-path[g,g]
  
 ind<-g %in% names(de)
 if (length(g)>0 & sum(ind)>=1) weight<-downstreamCpp(set, g, g[ind]) +1 else weight<-setNames(rep(0, length(g)),g)
wei<-setNames(rep(0, length(g)),g)
 wei[ind]<-weight 
 return(wei)
 }
 
 
 
prs<-function(all, de, pathways, nperm){

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