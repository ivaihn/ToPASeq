PCIscores <-function(expr, set, both.directions){

N<-length(nodes(set))
if (both.directions) g <- buildGraphNEL(nodes(set), edges(set), TRUE) else
          g <- buildGraphNEL(nodes(set), edges(set), FALSE)
a<-as(g,"matrix")
diag(a)<-1
a<-a[rownames(a) %in% rownames(expr), rownames(a) %in% rownames(expr), drop=FALSE]


if (nrow(which(a==1, TRUE))==0) {return(rep(NA, ncol(expr)))} else{
sc<-apply(which(a==1, TRUE),1, function(x) {
i<-x[1]
j<-x[2]
s<-sign(expr[i,]+ expr[j,])* sqrt(abs(expr[i,])) * sqrt(abs(expr[j,]))  
return(s)
})
if (length(dim(sc))==0) score=rep(0,ncol(expr)) else score<-rowSums(sc, na.rm=TRUE)
return(score/nrow(a))
}
}
