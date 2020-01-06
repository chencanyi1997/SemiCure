#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Multiply a number by two
//'
//' @param x A single integer.
//' @export
// [[Rcpp::export]]
double Inv_F_cpp(double u, double t0){
  double a = (1-u)*(1-exp(-t0));
  return -log(exp(-t0) + a);
}

// [[Rcpp::export]]
double G_gamma_cpp(double x, double gamma){
  return exp(-x)*(gamma==0) + pow((1+gamma*x),(-1/gamma))*(gamma>0);
}

// [[Rcpp::export]]
double Sfcn_cpp(double gamma, double F_hat, arma::vec beta, arma::vec Z){
  return G_gamma_cpp(exp(arma::dot(beta,Z))*F_hat,gamma);
}

