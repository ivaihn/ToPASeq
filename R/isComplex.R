isComplex<-function (graph, x) 
{   
    
    e <- edges(graph)
    sube <- e[e[, 1] %in% x & e[, 2] %in% x, ]
    nrow(sube) == ncol(combn(x,2)) & all(x %in% unlist(sube[, 1:2])) & 
        all(sube[, 3] == "undirected") & 
        all( (regexpr("binding",sube[, 4]) != -1) )
}

#example
#g<-convertIdentifiers(kegg[["p53 signaling pathway"]], "symbol")
#x<-c("CCNE1", "CDK2",  "CCNE2")
#isComplex(g,x)
#
#
#x<-c(x, nodes(g)[15])
#isComplex(g,x)

###########




