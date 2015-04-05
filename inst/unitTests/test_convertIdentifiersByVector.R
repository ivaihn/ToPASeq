test_convertIdentifiersByVector<-function(){
g<-kegg[["Asthma"]]
conv<-setNames(paste("gene", 1:length(nodes(g)), sep=""), nodes(g))
gc<-convertIdentifiersByVector(g, conv, "dummy")
checkTrue(length(nodes(g))== length(nodes(gc)))
checkTrue(gc@ident=="dummy")

conv<-conv[-1]
obs <- tryCatch(convertIdentifiersByVector(g, conv, "dummy"), warning=conditionMessage)
expe<-paste("These pathway nodes are missing in the 'conv.table':", paste(nodes(g)[1], collapse=", "),". The original identifiers were kept.")
checkIdentical(expe, obs)
}