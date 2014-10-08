processDEGraph<-function(res, overall){
if (length(res)>0) {
   nc<-max(unlist(sapply(res, length)))
   p<-suppressWarnings(t(sapply(res, function(x) sapply(1:nc, function(i) as.numeric(try(x[[i]][[1]][1], TRUE))))))
   pFourier<-suppressWarnings(t(sapply(res, function(x) sapply(1:nc, function(i) as.numeric(try(x[[i]][[1]][2], TRUE))))))
   graphs<-suppressWarnings(t(sapply(res, function(x) sapply(1:nc, function(i) 
   if ( any("graphNEL" == is((try(x[[i]][[2]], TRUE))))) x[[i]][[2]] else NA
 ))) )
   k<-suppressWarnings(t(sapply(res, function(x) sapply(1:nc, function(i) as.numeric(try(x[[i]][[3]], TRUE))))))

   ord<-as.numeric(sapply(1:nc, function(i) seq(i,4*nc, nc)))

out<-cbind(p, pFourier, graphs, k)[,ord]
tmp<-c("p","pFourier", "graph", "k")

colnames(out)<-paste("Comp",rep(1:nc, each=4),".",tmp, sep="")

if (!(any(overall==c("min","mean","biggest")))) stop("Invalid value for 'overall'. Use one of the following: 'min','mean','biggest'") 
if (overall=="min") f<-min 
if (overall=="mean") f<-mean 
if (overall=="biggest") f<-function(x, na.rm=TRUE) {x[1]}
overalp<-apply(pFourier, 1, function(x) f(x,na.rm=TRUE)) 

out<-cbind(Overall.p=overalp, out)  
}  else out<-res
return(out)
}
