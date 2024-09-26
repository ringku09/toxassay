// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP matrixProd(const Eigen::Map<Eigen::MatrixXd> A,
                Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
SEXP get_pi(const Eigen::Map<Eigen::MatrixXd> X,
            Eigen::Map<Eigen::MatrixXd> Y) {
  int m = X.rows();
  Eigen::VectorXd p(m);
  Eigen::MatrixXd V = X * Y;
  for (int i = 0; i < m; ++i) {
    p[i] = V.row(i).dot(X.row(i));
  }
  return Rcpp::wrap(p);
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
SEXP get_psi(const Eigen::Map<Eigen::MatrixXd> X,
             Eigen::Map<Eigen::SparseMatrix<double>> Y) {
  int m = X.rows();
  Eigen::VectorXd p(m);
  Eigen::MatrixXd V = X * Y;
  for (int i = 0; i < m; ++i) {
    p[i] = V.row(i).dot(X.row(i));
  }
  return Rcpp::wrap(p);
}


// SEXP hlm_cpp2(const Eigen::Map<Eigen::MatrixXd> X,
//              Eigen::Map<Eigen::MatrixXd> Y) {
//   Eigen::VectorXd p = ((X * Y).array() * X.array()).rowwise().sum();
//   return Rcpp::wrap(p);
// }
