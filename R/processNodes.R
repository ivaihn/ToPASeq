#processNodes <-
#function(set, de, all){
#if (!all(names(de) %in% all)) stop("Differentially expressed genes must be a subset of All genes")
#nodes.info<-sapply(nodes(set), function(x) {
# expr=0
# weight=0
# if (x %in% all)  expr=1
# if (x %in% names(de)) {expr=de[x]
# down<-downstream.nodes(x,set,names(de))
# weight=1+ sum(down %in% names(de))}
# return(c(expr, weight))
# })
#rownames(nodes.info)<-c("Log Fold-Change","Number of downstream DEG")
#return(nodes.info)
#}

# ifelse
processNodes <-function(set, de, all, down.list){
if (!all(names(de) %in% all)) stop("Differentially expressed genes must be a subset of All genes")

expr<-ifelse(nodes(set) %in% names(de), de[nodes(set)], ifelse(nodes(set) %in% all, 1, 0 ))
down.genes<-sapply(down.list[[set@title]], function(x) sum(x %in% names(de)))
weight<-ifelse(nodes(set) %in% names(de), 1+ down.genes ,0)

nodes.info<-rbind(expr, weight)
colnames(nodes.info)<-nodes(set)
rownames(nodes.info)<-c("Log Fold-Change","Number of downstream DEG")
return(nodes.info)
}
