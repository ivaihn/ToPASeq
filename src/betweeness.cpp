#include <Rcpp.h>
#include <queue>
#include <stack>
#include <iostream>



using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

int find(String s, CharacterVector CV) {
  int out=-1;            
  for (int i=0; i<CV.size(); i++) {
    String tmp=CV[i];
    if (tmp==s) {out=i; break;}
  }
  return(out);
}
// [[Rcpp::export]]

NumericVector betweenness(IntegerMatrix g, CharacterVector V) {
  
  int n=V.size();
  V.names()=V;
  NumericVector CB(n);
  CB.names()=V;
  for(int i=0; i<n; i++) {
    String s=V[i];
    stack <String> S;
    
    List P(n);
    P.names()=V;
    CharacterVector tmp;
    for (int ii=0; ii<n; ii++) P[ii]=tmp;
    
    NumericVector sigma(n);
    sigma.names() = V;
    sigma.fill(0);
    sigma[s]=1;
    
    NumericVector d(n);
    d.fill(-1);
    d.names() = V;
    d[s] = 0;
    
    queue <String> Q;
    Q.push(s);
    

    while (Q.size() >0 ){
      String v;
      v=Q.front();
      Q.pop();
      S.push(v);
      int iv =find(v,V) ;
      LogicalVector nei(V.size());
      nei = g(iv,_);
      CharacterVector neighbours=V[nei];
        for(int j=0; j<neighbours.size(); j++){
         String w=neighbours[j];
         if (d[w]<0) {
           Q.push(w);
           d[w] = d[v] + 1.0;
         }
         
         if (d[w]==d[v]+1.0) {
           sigma[w] = sigma[w] + sigma[v];
             CharacterVector tmp=P[w];
             tmp.push_back(v);
             P[w]=tmp;
         }
        }
        
     }
   NumericVector delta(n);
   delta.names()=V;
   while (S.size()>0) {
     String w = S.top();
     S.pop();
     CharacterVector Pw=P[w];
    
     for(int iii=0; iii<Pw.size(); iii++){
       String v=Pw[iii];
        delta[v]=delta[v]+sigma[v]/sigma[w]*(1.0+delta[w]);
      }
      if (w!=s) CB[w]=CB[w]+delta[w];
   }

  }
  return(CB);
}

