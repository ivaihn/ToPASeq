test_AdjacencyMatrix2pathway<-function(){

set.seed(1234)
adjmat<-matrix(rbinom(size=1, n=25, prob=0.3), 5,5)
rownames(adjmat)<-letters[1:5]
AM<-AdjacencyMatrix2pathway(adjmat)
checkEquals(nodes(AM), rownames(adjmat))
checkEquals(nlevels(edges(AM)[,4]),1)



colnames(adjmat)<-letters[6:10]
obs <- tryCatch(AdjacencyMatrix2pathway(adjmat), error=conditionMessage)
checkIdentical("Rownames and colnames of the matrix do not match", obs)

adjmat<-matrix(sample(c(0,0,0,0,1,-2),25, replace=TRUE ), 5,5)
colnames(adjmat)<-letters[6:10]
obs <- tryCatch(AdjacencyMatrix2pathway(adjmat), error=conditionMessage)
checkIdentical("The matrix contains values other than 0, 1, -1", obs)


adjmat<-matrix(sample(c(0,0,0,0,1,-1),25, replace=TRUE ), 5,5)
checkException(AdjacencyMatrix2pathway(adjmat))
obs <- tryCatch(AdjacencyMatrix2pathway(adjmat), error=conditionMessage)
checkIdentical("Adjacency matrix must have rownames and/or colnames", obs)

rownames(adjmat)<-letters[1:5]
AM<-AdjacencyMatrix2pathway(adjmat)
checkEquals(nodes(AM), rownames(adjmat))
checkEquals(nlevels(edges(AM)[,4]),2)

}