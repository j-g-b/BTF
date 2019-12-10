// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <vector>
#include <random>
#include <typeinfo>
#include <algorithm>
#include <iostream>
#include <array>        // std::array
#include <chrono>
#include <math.h>
#include <set>
#include <ctime>
#include <fstream>
//
using namespace arma;
using namespace Rcpp;
//
Eigen::MatrixXd RcppToEigenMatx(Rcpp::NumericMatrix & RMat){
  Eigen::MatrixXd EigMat = Eigen::MatrixXd::Zero(RMat.nrow(), RMat.ncol());
  for(int i = 0; i < RMat.nrow(); i = i + 1) {
    for(int j = 0; j < RMat.ncol(); j = j + 1){
      EigMat(i, j) = RMat(i, j);
    }
  }
  return(EigMat);
}
//
Eigen::VectorXd RcppToEigenVecx(Rcpp::NumericVector & RVec){
  Eigen::VectorXd EigVec = Eigen::VectorXd::Zero(RVec.size());
  for(int i = 0; i < RVec.size(); i = i + 1) {
    EigVec[i] = RVec[i];
  }
  return(EigVec);
}
//
Eigen::MatrixXd removeRow(Eigen::MatrixXd matrix, unsigned int j)
{
    unsigned int numRows = matrix.rows()-1;
    Eigen::MatrixXd new_mat = Eigen::MatrixXd::Zero(numRows, matrix.cols());
    if(j == 0){
        new_mat = matrix.bottomRows(numRows);
    } else if(j == 0){
        new_mat = matrix.topRows(numRows);
    } else {
        new_mat.topRows(j) = matrix.topRows(j);
        new_mat.bottomRows(numRows - j) = matrix.bottomRows(numRows - j);
    }
    //
    return new_mat;
}
//
Eigen::VectorXd removeVecRow(Eigen::VectorXd vector, unsigned int j)
{
    unsigned int numRows = vector.size()-1;
    Eigen::VectorXd new_vec = Eigen::VectorXd::Zero(numRows);
    if(j == 0){
        new_vec = vector.tail(numRows);
    } else if(j == 0){
        new_vec = vector.head(numRows);
    } else {
        new_vec.head(j) = vector.head(j);
        new_vec.tail(numRows - j) = vector.tail(numRows - j);
    }
    //
    return new_vec;
}
//
Eigen::VectorXd computeHatValues(Eigen::MatrixXd & XtXinv, Eigen::MatrixXd & X)
{
    Eigen::VectorXd hat_values = Eigen::VectorXd::Zero(X.rows());
    for(int i = 0; i < X.rows(); i = i + 1){
      hat_values[i] = X.row(i).dot(XtXinv*X.row(i).transpose());
    }
    //
    return hat_values;
}
//
// [[Rcpp::export]]
Rcpp::List rcpp_fabp_lin_reg(Rcpp::NumericVector & Y0, double & S0, Rcpp::NumericVector & R0, 
                             Rcpp::NumericMatrix & U0, Rcpp::NumericMatrix & V0)
{
  //
  // Set up Eigen objects
  Eigen::VectorXd Y = RcppToEigenVecx(Y0);
  Eigen::VectorXd R = RcppToEigenVecx(R0);
  Eigen::MatrixXd U = RcppToEigenMatx(U0);
  Eigen::MatrixXd V = RcppToEigenMatx(V0);
  //
  int N = U.rows();
  int M = V.rows();
  int D1 = U.cols();
  int D2 = V.cols();
  //
  double sigmasq_hat = S0;
  //
  // Form Kronecker product
  Eigen::MatrixXd X = Eigen::MatrixXd::Ones(M*N, D2*D1 + 1);
  for(int i = 0; i < M; i = i + 1){
    for(int j = 0; j < D2; j = j + 1){
        X.block(i*N, j*D1 + 1, N, D1) = V(i, j)*U;
    }
  }
  //
  Eigen::MatrixXd XtXinv = (X.transpose()*X).llt().solve(Eigen::MatrixXd::Identity(D2*D1 + 1, D2*D1 + 1));
  Eigen::MatrixXd XtXinvXt = XtXinv*X.transpose();
  //
  // Compute Prasad Rao estimator
  Eigen::VectorXd hatvals = computeHatValues(XtXinv, X);
  double pr_psi = ((Y - X*(XtXinvXt*Y)).array().square().matrix().sum() - (sigmasq_hat*(R.array().inverse())*(1 - hatvals.array())).matrix().sum()) / (M*N - D1*D2 - 1);
  //
  // Compute linking model estimators
  Eigen::VectorXd W = (sigmasq_hat*(R.array().inverse()) + pr_psi).inverse().matrix();
  Eigen::VectorXd Theta_hat = Eigen::VectorXd::Zero(M*N);
  Eigen::VectorXd Theta_tild = Eigen::VectorXd::Zero(M*N);
  Eigen::VectorXd Tau_tild = Eigen::VectorXd::Zero(M*N);
  //
  XtXinv = (X.transpose()*W.asDiagonal()*X).llt().solve(Eigen::MatrixXd::Identity(D2*D1 + 1, D2*D1 + 1));
  hatvals = computeHatValues(XtXinv, X);
  XtXinvXt = XtXinv*X.transpose();
  Theta_hat = X*XtXinvXt*W.asDiagonal()*Y;
  //
  for(int j = 0; j < N*M; j = j + 1){
    //
    Theta_tild[j] = Theta_hat[j] - (W[j]*hatvals[j]*(Y[j] - Theta_hat[j]))/(1 - W[j]*hatvals[j]);
    Tau_tild[j] = sigmasq_hat*(hatvals[j] / (1 - W[j]*hatvals[j])) + pr_psi;
  }
  //
  Rcpp::NumericVector rTau_tild(Tau_tild.size());
  Rcpp::NumericVector rTheta_tild(Theta_tild.size());
  //
  for(int i = 0; i < Tau_tild.size(); i = i + 1){
      rTau_tild[i] = Tau_tild[i];
  }
  //
  for(int i = 0; i < Theta_tild.size(); i = i + 1){
      rTheta_tild[i] = Theta_tild[i];
  }
  //
  return Rcpp::List::create(Rcpp::Named("Tau_tild") = rTau_tild,
                            Rcpp::Named("Theta_tild") = rTheta_tild);
  //
}
