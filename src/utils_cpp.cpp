#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// // [[Rcpp::export]]
// double Inv_F_generic(double u, double t0){
//   double a = (1-u)*(1-exp(-t0));
//   return -log(exp(-t0) + a);
// }

// [[Rcpp::export]]
double G_gamma_cpp(double x, double gamma){
  return exp(-x)*(gamma==0) + pow((1+gamma*x),(-1/gamma))*(gamma>0);
}

// // [[Rcpp::export]]
// double Sfcn(double gamma, double F_hat, arma::vec beta, arma::vec Z){
//   return G_gamma(exp(arma::dot(beta,Z))*F_hat,gamma);
// }

// [[Rcpp::export]]
double G_inv_cpp(double y, double gamma){
  return (-log(y)*(gamma == 0) + pow(y,(-gamma) - 1)/gamma*(gamma>0));
}
//
// // // [[Rcpp::export]]
// // arma::vec F_hat(arma::vec alpha){
// //   int m = alpha.n_elem + 1;
// //   arma::vec res = arma::zeros(m+1,1);
// //   double sum = 0;
// //   for(int k = 1; k < (m+1); k++){
// //     sum = sum + exp(alpha(k-1));
// //     res.row(k) = 1 - exp(-sum);
// //   }
// //   res.row(0) = 0;
// //   res.row(m) = 1;
// //   return res;
// // }
//
