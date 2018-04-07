#include <Rcpp.h> 
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector colMax(NumericMatrix X) {
    int ncol = X.ncol();
    Rcpp::NumericVector out(ncol);
    for (int col = 0; col < ncol; col++){
        out[col]=Rcpp::max(X(_, col)); 
    } 
    return wrap(out);
}

// [[Rcpp::export]]
NumericVector colMin(NumericMatrix X) {
    int ncol = X.ncol();
    Rcpp::NumericVector out(ncol);
    for (int col = 0; col < ncol; col++){
        out[col]=Rcpp::min(X(_, col)); 
    } 
    return wrap(out);
}
