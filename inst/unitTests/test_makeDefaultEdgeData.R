test_makeDefaultEdgeData<-function(){
edgdata<-makeDefaultEdgeData()
checkEquals(length(edgdata),2)
checkEquals(names(edgdata),c("graphite2SPIA","beta"))
checkEquals(colnames(edgdata[[1]]),c("type","spiaType"))
checkTrue(all(edgdata[[1]][,2] %in% edgdata[[2]][,1]))
checkTrue(all(unique(edgdata[[2]][,2]) %in% c(-1,1,0)))
}