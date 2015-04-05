#x draha triedy pathway


#x GEDM, rows=genes, col=samples
normalizationTAPPA<-function(x){
if (is.null(dim(x))) stop("Only matrices can be normalized")
x<-apply(x, 2, function(y) (y-mean(y))/sd(y)) #colum to zero mean and same scope
x<- 1/(1+exp(-x))-0.5 
return(x)
}

#x normalized GEDM, rows=genes col=samples
#A adjacency matrix
calculatePCI<-function(x, A){
sgn<-function(x){ifelse(x>0, 1, ifelse(x<0, -1, 0))}
# je normalizovana matica?
# je A spravna matica - s 1 na diagonale

#
A.ind<-which(A==1, arr.ind=TRUE)
if (nrow(A.ind>0)) {
PCI<-apply(A.ind, 1, function(ind) {
i<-ind[1]
j<-ind[2]
PCI<-sgn(x[i,]+x[j,])*sqrt(abs(x[i,]))*A[i,j]*sqrt(abs(x[j,]))
})
PCI<-rowSums(PCI)
norm.PCI=PCI/nrow(A)
} else norm.PCI<-NA
return(norm.PCI)
}


extractsubset<-function(x,A){
if (is.matrix(A)) {
valid.genes<-rownames(A) %in% rownames(x)
if (sum(valid.genes)==0) warning("Expression data do not match node labels")
return(list(x=x[rownames(A)[valid.genes],], A=A[valid.genes, valid.genes]))
}
if (any(is(A)=="graphNEL")) {
valid.genes<-nodes(A) %in% rownames(x)
if (sum(valid.genes)==0) warning("Expression data do not match node labels")
return(list(x=x[nodes(A)[valid.genes],], A=A))
}

}


TAPPASingle<-function(path, group,x, test,  verbose) {
if (verbose) message(path[[2]])
A<-path[[1]]
valid.data<-extractsubset(x,A)
PCI<-calculatePCI(valid.data$x,valid.data$A)
desc<-unlist(tapply(PCI, group, function(x) c(N=length(x[!is.na(x)]), summary(x, na.rm=TRUE))))
out<-c(desc, p.value=test(PCI~group)$p.value)
return(out)
}


tappa<-function(x, group, pathways, test, normalize=TRUE, verbose=FALSE ){

if (normalize) x<-normalizationTAPPA(x)

out<-catchErr(pathways, function(p) TAPPASingle(p, group, x,test, verbose))

if (length(out[[1]])>0) {
out[[1]]<-data.frame(t(vapply(out[[1]], function(x) x, numeric(15))))
out[[1]]$q.value<-p.adjust(out[[1]]$p.value,"fdr")
}
return(out)
}



