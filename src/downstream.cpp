#include <Rcpp.h>

using namespace Rcpp;
// [[Rcpp::export]]

IntegerVector downstreamCpp(    IntegerMatrix AM, 
                                CharacterVector nodeNames, 
                                CharacterVector index   )
{
int nN = nodeNames.size();
int nIndex = index.size();
IntegerVector whN = match(index, nodeNames);

bool miss = is_true(any(is_na(whN)));

if (miss) stop ("unmatched node provided");

 IntegerVector rval(index.size(),0);  // len celkovy pocet bez konkretnych vysledkov
 for (int i=0; i < nIndex; i++) {
 
  IntegerVector marked(nN);
  IntegerVector distv(nN);
  marked.names() = nodeNames;
  distv.names() = nodeNames;
  int distx = 1;
  int newmk = 0;
  String current = index[i];
  marked[current] = 1;
  bool done = false;
  while(!done) {
  
    LogicalVector selN = marked==1;
    CharacterVector minds = nodeNames[selN];
    //Rf_PrintValue(marked);
    //Rf_PrintValue(minds);
    for (int j=0; j<minds.size(); j++ ) {
    // CharacterVector  nodeS = minds(j);
    int node = match(minds, nodeNames)(j);
    IntegerVector avec = AM(node-1,_);
    avec.names() = nodeNames;
    LogicalVector discard = is_na(match(nodeNames,index));
    avec[discard] = 0;
    LogicalVector discard2 = (marked != 0) & (avec == 1);
    avec[discard2] = 0;
    LogicalVector avecFin = avec==1;
    CharacterVector avecNames = nodeNames[avecFin];
    //Rf_PrintValue(avec);
    marked[avecNames] = 1;
    distv[avecNames] = distx;
    }
    
    marked[minds] = 2;
    distx++;
    newmk = sum(marked == 1);
    if (newmk == 0) done = true;
    }
  marked[current] = 0;
  rval[i] = sum(marked==2);
  }
  

 rval.names()=index;
 return(rval);
 }



