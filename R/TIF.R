TIF <-
function(expr, set, both.directions, alpha=0.05){
g<-rownames(expr)
if (both.directions) set <- buildGraphNEL(nodes(set), edges(set), TRUE) else
          set <- buildGraphNEL(nodes(set), edges(set), FALSE)
          
sp<-shortest.paths(igraph.from.graphNEL(set))
select.nodes<-nodes(set)%in% g
sp<-sp[select.nodes, select.nodes] 

pc<-cor(t(expr), use="pairwise.complete.obs")
psi<-exp(-sp/abs(pc))

#n<-length(nodes(set))
f<-sp/abs(pc)
diag(f)<-0


tif<-sapply(1:nrow(f), function(i) {
    x<-f[i,]
    x<-x[f[i,] >= -log(alpha)]
    x<-x[is.finite(x)]
    if (length(x)==0) return(1) else return(1+exp(-mean(x )))
})
names(tif)<-rownames(expr)
return(tif)
}
