clipperSingle<-function(expr, classes, g, nperm, alphaV, b, permute, method, easyClip){

if (easyClip) 
  out<-list(unlist(pathQ(expr, classes, g, nperm=nperm, alphaV=alphaV, b=nperm, permute=TRUE)),
          easyClip(expr, classes, g, method =method,nperm=nperm, b=nperm)) else
  out<-list(unlist(pathQ(expr, classes, g, nperm=nperm, alphaV=alphaV, b=nperm, permute=TRUE)),
            list())

 return(out)
}  
        

CLIPPER<-function (pathways, expr, classes, method, testCliques=FALSE, nperms, alphaV=0.05, b=NULL, permute=TRUE, easyClip=FALSE){
    #if (!require(clipper))        stop("library clipper is missing")

    classes<-factor(as.numeric(factor(classes)))
message("Analysing pathway:\n")
      out<-catchErr(pathways, function(p) {
      cat(p[[2]],"\n")
clipperSingle(expr, classes, p[[1]], nperms, alphaV, b, permute, method, easyClip)
})
cliq.test<-list()
if (testCliques) {
message("Testing cliques\n")
 cliq.test<-lapply(pathways, function(p) {
 cat(p[[2]],"\n")
 if (method=="mean") cliq<-cliqueMeanTest(expr, classes, p[[1]], nperms, alphaV, b, permute=permute)
 if (method=="var") cliq<-cliqueVarianceTest(expr, classes, p[[1]], nperms, alphaV, b, permute=permute)
 return(cliq)
 })
 
  }
if (length(out[[1]])>0) {
res<-data.frame(t(sapply(out[[1]],function(x) x[[1]])))
res$mean.q.value<-p.adjust(res$alphaMean,"fdr")
res$var.q.value<-p.adjust(res$alphaVar,"fdr")

paths<-lapply(out[[1]],function(x) x[[2]])

out[[1]]<-list(results=res, paths=paths, cliq.test)
} else 
out[[1]]<-list(results=out[1:2], paths=list(), cliq.test)
return(out)
}

CLIPPERfast<-function (pathways, expr, classes, nperms, alphaV=0.05, b=NULL, permute=TRUE){
    #if (!require(clipper))        stop("library clipper is missing")

    clipperSingle<-function(expr, classes, g, nperm, alphaV){

out<-list(unlist(pathQ(expr, classes, g, nperm=nperm, alphaV=alphaV, b=nperm, permute=TRUE)))
 return(out)
} 
    classes<-factor(as.numeric(factor(classes)))
message("Analysing pathway:\n")
      out<-catchErr(pathways, function(p) {
      cat(p[[2]],"\n")
clipperSingle(expr, classes, p[[1]], nperms, alphaV)
})
cliq.test<-list()

if (length(out[[1]])>0) {
res<-data.frame(t(sapply(out[[1]],function(x) x[[1]])))
res$mean.q.value<-p.adjust(res$alphaMean,"fdr")
res$var.q.value<-p.adjust(res$alphaVar,"fdr")

out[[1]]<-list(res)
}
return(out)
}
    
