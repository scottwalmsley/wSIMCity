#include <Rcpp.h>
using namespace Rcpp;

//' @name ppmErr
//' @title Calculate ppm error for a given mass
//'
//' @description This function compares 2 masses and computes the ppm error between them,
//' @param m a double number, the experimental mass
//' @param m0 a double number, the theoretical mass
//' @export
//' @return double number the ppm error
//' @examples
//' ppmErr(490.1944,490.1946)
//'  -0.4080012
// [[Rcpp::export]]
double ppmErr(double m, double m0) {
  return (m-m0)/m0*1e6;
}


//' @name getMassTolRange
//' @title Get mass tolerance range in ppm
//'
//' @description  Generates the tolerance window in Da for a given ppm value and mass
//' @param m a double number, the mass
//' @param ppm a double number, the user provided ppm tol window
//' @export
//' @return vector of numeric values containing min and max range in Da
//' @examples
//' getMassTolRange(m = 490.1946,ppm = 5)
//' 490.1921 490.1971
// [[Rcpp::export]]
NumericVector getMassTolRange(long double m,long double ppm){
  double dM = ppm*m/1e6;
  NumericVector out = NumericVector::create(m-dM,m+dM);
  
  return out;
}


//' @name dlaplace
//' @title Probability density function: Laplace
//'
//' @description  Generates the tolerance window in Da for a given ppm value and mass
//' @param X numeric vector of values
//' @param m numeric the mode of the distribution
//' @param b numeric scale parameter
//' @export
//' @return vector of numeric values containing density values
//' @examples
//' x <- seq(-10,10,by = 0.001)
//' d_lap <- dlaplace(x,m=0,b = 1)
//' plot(x,d_lap, type = "l")
// [[Rcpp::export]]
NumericVector dlaplace( NumericVector X, double m, double b) {
  
  NumericVector out(X.size());
  
  for(int i; i<X.size();i++){
   
      out[i] =   1.0 / (2.0 * b) *std::exp(-1*std::abs(X[i] - m) / b);

  }
  
  return out;

}

//' @name weight_laplace
//' @title Calculate weight with Laplace assumption.
//'
//' @description  Assigned a weight to the distribution of X values given expected instrument performance.
//' For example, dM is the input X values, tol is a tolderance in ppm, and  boost is the multiplier.
//' @param dM vector or numeric the vector of measurement errors
//' @param tol numeric the ppm tolerance of the instrument
//' @param boost numeric the boost value (0-2)
//' @param mu numeric the mean of the laplace distribution
//'
//' @return numeric or vector the weighted value(S).
//' @export
//'
//' @examples 
//' weight_laplace(dM,mu=0.1, tol = 5, boost = 2)
// [[Rcpp::export]]
NumericVector weight_laplace(NumericVector dM, double mu,double tol,double boost){
  
  return boost*exp(-1* (abs(mu-dM)/tol));

}




