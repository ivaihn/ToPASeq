normalizeTAPPA <-
function(dat){
dat<-apply(dat,2, function(x) {
 x<-(x-mean(x))/sd(x) 
 x<-1/(1+exp(-x))-0.5 
})
return(dat)
}
