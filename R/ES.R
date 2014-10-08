ES <-
function(score1,score2){
all<-c(score1, score2)
rank.all<-rank(-all)
ranks.in<-rank.all[1:length(score1)]
ranks.out<-rank.all[(length(score1)+1):length(all)]


sc1<-rep(0,length(all))
sc2<-rep(0,length(all))
sc1[ranks.in]<-sort(score1)/sum(score1)
n.out=length(score2)
sc2[ranks.out]<-1/n.out
#score1<-sort(score1)
P.in=cumsum(sc1)
P.out=cumsum(sc2)
es<-max(abs(range(P.in-P.out)))#[which.max(abs(range(P.in-P.out)))]
return(es)
}
