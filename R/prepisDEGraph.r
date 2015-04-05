
processDEGraph<-function(res, overall){
if (length(res)>0) {
   nc<-max(unlist(sapply(res, length)))
   p<-suppressWarnings(t(sapply(res, function(x) sapply(1:nc, function(i) as.numeric(try(x[[i]][[1]][1], TRUE))))))
   pFourier<-suppressWarnings(t(sapply(res, function(x) sapply(1:nc, function(i) as.numeric(try(x[[i]][[1]][2], TRUE))))))
   graphs<-suppressWarnings(t(sapply(res, function(x) sapply(1:nc, function(i) 
   if ( any("graphNEL" == is((try(x[[i]][[2]], TRUE))))) x[[i]][[2]] else NA
 ))) )
   colnames(graphs)<-paste("Comp",1:nc,".","graph", sep="")
   k<-suppressWarnings(t(sapply(res, function(x) sapply(1:nc, function(i) as.numeric(try(x[[i]][[3]], TRUE))))))

   ord<-as.numeric(sapply(1:nc, function(i) seq(i,3*nc, nc)))

out<-cbind(p, pFourier, k, deparse.level=0)[,ord]
tmp<-c("p","pFourier", "k")

colnames(out)<-paste("Comp",rep(1:nc, each=3),".",tmp, sep="")



if (!(any(overall==c("min","mean","biggest")))) stop("Invalid value for 'overall'. Use one of the following: 'min','mean','biggest'") 
if (overall=="min") f<-min 
if (overall=="mean") f<-mean 
if (overall=="biggest") f<-function(x, na.rm=TRUE) {x[1]}
overalp<-apply(pFourier, 1, function(x) f(x,na.rm=TRUE)) 

out<-cbind(Overall.p=overalp, out)  
out<-list(data.frame(out), graphs=graphs)
}  else out<-res
return(out)
}

degraph<-function(exprs, group,pathways,  overall="biggest"){

if (!require(DEGraph)) {
        stop("library DEGraph is missing")
    checkPkgVersion("DEGraph", "1.4.0")
   } else {
out<-catchErr(pathways, function(p) {
if (require(DEGraph)) DEGraph::testOneGraph(p[[1]], exprs, group, useInteractionSigns = FALSE)
})

out[[1]]<-processDEGraph(out[[1]], overall)
out[[1]][[1]]$Overall.q.value<-p.adjust(out[[1]][[1]][,"Overall.p"],"fdr")
out[[1]][[1]]<-out[[1]][[1]][,c(1, ncol(out[[1]][[1]]), 2:(ncol(out[[1]][[1]])-1))]
return(out)
}}