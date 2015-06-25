getdatp<-function(x, rel, Beta){
names(Beta)<-rel
con<-Reduce("+",Map(function(m, b) {m*abs(sign(b))}, x[rel],as.list(Beta[rel])), init=0  )
s<- Reduce("+",Map(function(m, b) {m*b}, x[rel],as.list(Beta[rel])), init=0  )
z = matrix(rep(apply(con, 2, sum), dim(con)[1]), dim(con)[1], dim(con)[1], byrow = TRUE)
z[z == 0] <- 1
return(s/z)
}

checkDEandAll<-function(de, all){
    IDsNotP <- names(de)[!names(de) %in% all]
    if (length(IDsNotP)/length(de) > 0.01) {
        stop("More than 1% of your de genes have IDs are not present in the reference array!. Are you sure you use the right reference array?")
    }
    if (!length(IDsNotP) == 0) {
        cat("The following IDs are missing from all vector...:\n")
        cat(paste(IDsNotP, collapse = ","))
    }
    if (length(intersect(names(de), all)) != length(de)) {
        stop("de must be a vector of log2 fold changes. The names of de should be included in the refference array!")
    }
}    
    
combfunc<-function (p1 = NULL, p2 = NULL, combine = "fisher") {
    tm = na.omit(c(p1, p2))
    if (!all(tm >= 0 & tm <= 1)) {
        stop("values of p1 and p2 have to be >=0 and <=1 or NAs")
    }
    if (combine == "fisher") {
        k = p1 * p2
        comb = k - k * log(k)
        comb[is.na(p1)] <- p2[is.na(p1)]
        comb[is.na(p2)] <- p1[is.na(p2)]
        return(comb)
    }
    if (combine == "norminv") {
        comb = pnorm((qnorm(p1) + qnorm(p2))/sqrt(2))
        comb[is.na(p1)] <- p2[is.na(p1)]
        comb[is.na(p2)] <- p1[is.na(p2)]
        return(comb)
    }
}

    
   
    ############
    # M jedna z datp matic
SPIASingle<-function(de, all, M, nB, combine="fisher"){

 ph <- pb <- pcomb <- nGP <- pSize <- smPFS <- tA <- tAraw <- NULL

diag(M) <- diag(M) - 1
X <- de[rownames(M)]
noMy <- sum(!is.na(X))
okg <- intersect(rownames(M), all)
ok <- rownames(M) %in% all

if (!((noMy) > 0 & (abs(det(M)) > 1e-07))) {
  pb <- ph <- smPFS <- pcomb <- tAraw <- tA <- ob <- NA 
  pfs<- rep(NA, length(X)) } else
{
X[is.na(X)] <- 0
pfs <- solve(M, -X)



pfstmp<-replicate(nB, {
x <- rep(0, length(X))
names(x) <- rownames(M)
x[ok][sample(1:sum(ok), noMy)] <- as.vector(sample(de,noMy))
tt <- solve(M, -x)
sum(tt-x)
})


mnn <- median(pfstmp)
pfstmp <- pfstmp - mnn
ob <- sum(pfs - X) - mnn


if (ob > 0) pb <- sum(pfstmp >= ob)/length(pfstmp) * 2
if (ob < 0) pb <- sum(pfstmp <= ob)/length(pfstmp) * 2
if (ob == 0) if (all(pfstmp == 0)) pb <- NA else pb <- 1

if (!is.na(pb)) {
if (pb <= 0)  pb <- 1/nB/100 
if (pb > 1)   pb <- 1 
                }
           
}


nGP   <- noMy
pSize <- length(okg)             
smPFS <- round(sum(pfs - X),3)
tAraw <- round(smPFS,3)
ph    <- round(phyper(q = noMy - 1, m = length(okg), n = length(all) - length(okg), k = length(de), lower.tail = FALSE),3)
tA    <- round(ob,3)
pcomb <- round(combfunc(pb, ph, combine),3)

return(c(NDE=nGP, pSize=pSize, smPFS=smPFS, tAraw=tAraw, pNDE=ph, tA=tA, pG=pcomb, pPERT=pb))
}

#pathways - datp matice, po kontrole na beta, rel
spia<-function(de, all, pathways, perm, combine){

checkDEandAll(de, all)   

out<-catchErr(pathways, function(p) SPIASingle(de, all, p[[1]], perm, combine))
out[[1]]<-data.frame(t(sapply(out[[1]], function(x) x)))

if (length(out[[1]])>0) {
tmp<-out[[1]]  
 
    tmp$pGFdr = p.adjust(tmp$pG, "fdr")
    tmp$pGFWER = p.adjust(tmp$pG, "bonferroni")
    tmp$Status = ifelse(tmp$tA > 0,"Activated", "Inhibited"  )

tmp<-tmp[,c("pSize", "NDE", "pNDE", "tA", "pPERT", "pG", "pGFdr", "pGFWER", "Status")]
out[[1]]<-tmp    
}
return(out)
}


SPIAweights<-function(de, all, M){

diag(M) <- diag(M) - 1
X <- de[rownames(M)]
noMy <- sum(!is.na(X))
okg <- intersect(rownames(M), all)
ok <- rownames(M) %in% all

if (!((noMy) > 0 & (abs(det(M)) > 1e-07))) {
 pfs<-setNames(rep(NA, length(X)), rownames(M))
 
 } else {
X[is.na(X)] <- 0
pfs <- solve(M, -X)
}
return(pfs-X)
}

collectWeightsSPIA<-function(de, all, pathways){
out<-catchErr(pathways, function(p) SPIAweights(de, all, p[[1]]))
return(out[[1]])
}



SPIA4plot<-function(de,all,M, nB, name, combine="fisher"){


 ph <- pb <- tA <-  NULL

diag(M) <- diag(M) - 1
X <- de[rownames(M)]
noMy <- sum(!is.na(X))
okg <- intersect(rownames(M), all)
ok <- rownames(M) %in% all

if (!(noMy) > 0 & (abs(det(M)) > 1e-07)) {
  pb <- ph  <- tA <- NA } else
{
X[is.na(X)] <- 0
pfs <- solve(M, -X)

pfstmp<-replicate(nB, {
x <- rep(0, length(X))
names(x) <- rownames(M)
x[ok][sample(1:sum(ok), noMy)] <- as.vector(sample(de,noMy))
tt <- solve(M, -x)
sum(tt-x)
})


mnn <- median(pfstmp)
pfstmp <- pfstmp - mnn
ob <- sum(pfs - X) - mnn


if (ob > 0) pb <- sum(pfstmp >= ob)/length(pfstmp) * 2
if (ob < 0) pb <- sum(pfstmp <= ob)/length(pfstmp) * 2
           
if (pb <= 0)  pb <- 1/nB/100 
if (pb > 1)   pb <- 1 
                
if (ob == 0 & all(pfstmp == 0)) pb <- NA  else pb <- 1
           
}

tA    <- ob 
return(list(log2FCH=X, Pertubation=pfs, NetPertubation=pfs-X, RandomPertubation=pfstmp, pPERT=pb, ObservedTotalPERT=tA))}




drawSPIA<-function(x, name){
X<-x$log2FCH
pfs<-x$Pertubation
pfstmp<-x$RandomPertubation
pb<-x$pPERT
tA<-x$ObservedTotalPERT
xr<-range(X)
yr<-range(pfs-X)
par(mfrow = c(1, 2))
plot(X, pfs - X, xlim=xr, ylim=yr, main = paste("pathway ID=", name, sep = ""), xlab = "Log2 FC", ylab = "Perturbation accumulation (Acc)", cex.main = 0.8, cex.lab = 1.2)
abline(h = 0, lwd = 2, col = "darkgrey")
abline(v = 0, lwd = 2, col = "darkgrey")
points(X[abs(X) > 0 & X == pfs], pfs[abs(X) > 0 & X == pfs] - X[abs(X) > 0 & X == pfs], col = "blue", pch = 19, cex = 1.4)
points(X[abs(X) > 0 & X != pfs], pfs[abs(X) > 0 & X != pfs] - X[abs(X) > 0 & X != pfs], col = "red", pch = 19, cex = 1.4)
points(X[abs(X) == 0 & X == pfs], pfs[abs(X) == 0 & X == pfs] - X[abs(X) == 0 & X == pfs], col = "black", pch = 19, cex = 1.4)
points(X[abs(X) == 0 & X != pfs], pfs[abs(X) == 0 & X != pfs] - X[abs(X) == 0 & X != pfs], col = "green", pch = 19, cex = 1.4)
legend(x="topright", legend=c("DE, not pertubed","DE, pertubed","nonDE, not pertubed","nonDE, pertubed"), col=c("blue", "red", "black","green"), pch=19)
if (require(plotrix)) plotrix::thigmophobe.labels(X, pfs-X, labels=names(X)) else text(X, pfs-X, labels=names(X), pos=3)
plot(density(pfstmp, bw = sd(pfstmp)/4), cex.lab = 1.2, col = "black", lwd = 2, main = paste("pathway ID=", name, "  P PERT=", 
   round(pb, 5), sep = ""), xlim = c(min(c(tA - 0.5, pfstmp)), max(c(tA + 0.5, pfstmp))), cex.main = 0.8, 
   xlab = "Total Perturbation Accumulation (TA)")
abline(v = 0, col = "grey", lwd = 2)
abline(v = tA, col = "red", lwd = 3)

}    
    ############
    
    

 