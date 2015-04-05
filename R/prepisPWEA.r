#TIF
#path = graphNEL draha
#x - expresia genov z drahy, rows=genes
TIF<-function(path, x, alpha){
d<- johnson.all.pairs.sp(path)
pcc<-cor(t(x), method="pearson", use="pairwise.complete.obs")
d<-d[rownames(pcc),rownames(pcc)]
f<-d/abs(pcc)
diag(f)<-NA
valid<- f <= -log(alpha) 
f<-as.data.frame(t(f))
valid<-as.data.frame(t(valid))
tif<-unlist(Map(function(x,y) if (sum(y, na.rm=TRUE)>0) mean(x[y], na.rm=TRUE) else NA, f,valid))
tif<-exp(-tif)
tif[is.na(tif)]<-0

return(tif)
}


#tif list tif vsetkych genov zo vsetkych drah 
# alpha - alpha z povodnej metody
# n - pocet genov mimo drahy
notPTIF<-function(tif, n){
tif.all<-unlist(tif)
prop.pass<-sum(tif.all>1)/length(tif)

tif.pass<-tif.all[tif.all > 0]

input<-rbinom(n,size=1,prob=prop.pass)
input[input==1]<-rnorm(sum(input), mean=mean(tif.pass), sd=sd(tif.pass))
return(input)
}

#geneList gene level statistics (t) ORDERED
# tif of all genes (exponent)
# g genes in pathway (geneSet)
#  zbalika DOSE, fortify = vsetky hodnoty (pre graf)
gseaScores<-function (geneList, geneSet, exponent = 1, fortify = FALSE) 
{
    geneSet <- intersect(geneSet, names(geneList))
    N <- length(geneList)
    Nh <- length(geneSet)
    Phit <- Pmiss <- numeric(N)
    hits <- names(geneList) %in% geneSet
    Phit[hits] <- abs(geneList[hits])^exponent
    NR <- sum(Phit)
    Phit <- cumsum(Phit/NR)
    Pmiss[!hits] <- 1/(N - Nh)
    Pmiss <- cumsum(Pmiss)
    runningES <- Phit - Pmiss
    max.ES <- max(runningES)
    min.ES <- min(runningES)
    if (abs(max.ES) > abs(min.ES)) {
        ES <- max.ES
    }    else {
        ES <- min.ES
    }
    if (fortify == TRUE) {
        df <- data.frame(x = seq_along(runningES), runningScore = runningES, 
            position = as.integer(hits))
        return(df)
    }
    return(ES)
}

gseaScoresCols<-function (geneList, geneSet, exponent = 1) {
    geneSet <- intersect(geneSet, rownames(geneList))
    N <- nrow(geneList)
    Nh <- length(geneSet)

    Phit <- Pmiss <- matrix(0, nrow=N, ncol=ncol(geneList))

    hits <- rownames(geneList) %in% geneSet

    Phit[hits,] <- abs(geneList[hits,])^exponent
    NR <- colSums(Phit)
    NR <- matrix(rep(NR, each=N), nrow=N)
    Phit <- apply(Phit/NR, 2, cumsum)
    Pmiss[!hits,] <- 1/(N - Nh)
    Pmiss <- apply(Pmiss, 2, cumsum)
    runningES <- Phit - Pmiss

    max.ES <- colMax(runningES)
    min.ES <- colMin(runningES)
    ES<-ifelse(abs(max.ES) > abs(min.ES), max.ES, min.ES)
   
    return(ES)
}

PWEAscore<-function(tif, tifout, geneList, geneSet){

genPow<-function(x, pow){ return(sign(x)*abs(x)^(pow)) }

if (length(geneList) != length(tif)+length(tifout)) stop("Number of Topology Impact Factors differs from the number of gene-level statistics")

tif.complete<-vector("numeric", length(geneList))
names(tif.complete)<-names(geneList)

tif.complete[names(tif)]<-tif
tif.complete[names(tifout)]<-tifout

r<-genPow(abs(geneList),(1+tif.complete))
r<-sort(r, decreasing=TRUE)

sc<-gseaScores(r, geneSet, 1+tif)
return(sc)
}

preparePermsPWEA<-function(x, gr, nperm, test){
out<-replicate(nperm, test(x, sample(gr))$stat)
return(out)
}

prepareTIF<-function(pathways, exprs, alpha){

all.genes<-rownames(exprs)

inP<-lapply(pathways, function(x) {
 exprs.valid<-extractsubset(exprs, x)$x
 return(TIF(x, exprs.valid, alpha))
 }
 )

outP<-lapply(pathways, function(p){
g<-sum(! all.genes %in% nodes(p))
return(notPTIF(inP, g))
})

return(Map(c, inP, outP))
}

PWEASingle<-function(p, geneList, tif, perms){
obs<-gseaScores(geneList, nodes(p), tif, fortify=FALSE)
rnd<-gseaScores(perms, nodes(p), tif)
p<-sum(obs>= rnd)/ length(rnd)
return(c(ES=obs, p.value=p))
}

pwea<-function(obs, perms, tif, pathways, alpha){

geneList<-obs$stat

out<-catchErr(pathways, function(p) PWEASingle(p[[1]], geneList, tif[p@title], perms))

out[[1]]$q.value<-p.adjust(out[[1]]$p.value,"fdr")
return(out)
}
