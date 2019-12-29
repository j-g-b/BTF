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
  // Permute rows of X, Y, and R
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(N*M);
  perm.setIdentity();
  std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size());
  X = perm * X;
  Y = perm * Y;
  R = perm * R;
  //
  Eigen::MatrixXd XtXinv;
  Eigen::MatrixXd XtXinvXt;
  Eigen::VectorXd hatvals;
  double pr_psi;
  Eigen::VectorXd Theta_hat;
  Eigen::VectorXd Theta_tild = Eigen::VectorXd::Zero(M*N);
  Eigen::VectorXd Tau_tild = Eigen::VectorXd::Zero(M*N);
  Eigen::VectorXd W;
  Eigen::MatrixXd Xblock;
  Eigen::VectorXd Psi_tild = Eigen::VectorXd::Zero(2);
  Eigen::VectorXd Yseg;
  Eigen::VectorXd Rseg;
  //
  Eigen::MatrixXd blocks(2, 2);
  blocks(0, 0) = 0;
  blocks(0, 1) = N*M/2 - 1;
  blocks(1, 0) = N*M/2;
  blocks(1, 1) = N*M;
  //
  for(int b = 0; b < 2; b = b + 1){
    //
    // Compute Prasad Rao estimator
    Xblock = X.block(blocks(b, 0), 0, N*M / 2, X.cols());
    XtXinv = (Xblock.transpose()*Xblock).llt().solve(Eigen::MatrixXd::Identity(X.cols(), X.cols()));
    XtXinvXt = XtXinv*Xblock.transpose();
    hatvals = computeHatValues(XtXinv, Xblock);
    Yseg = Y.segment(blocks(b, 0), N*M / 2);
    Rseg = R.segment(blocks(b, 0), N*M / 2);
    pr_psi = ((Yseg - Xblock*(XtXinvXt*Yseg)).array().square().matrix().sum() - (sigmasq_hat*(Rseg.array().inverse())*(1 - hatvals.array())).matrix().sum()) / ((M*N/2) - D1*D2 - 1);
    pr_psi = std::max(pr_psi, 0.0);
    //
    // Compute linking model estimators
    Theta_hat = Eigen::VectorXd::Zero(M*N);
    //
    W = (sigmasq_hat*(R.array().inverse()) + pr_psi).inverse().matrix();
    XtXinv = (X.transpose()*W.asDiagonal()*X).llt().solve(Eigen::MatrixXd::Identity(X.cols(), X.cols()));
    hatvals = computeHatValues(XtXinv, X);
    XtXinvXt = XtXinv*X.transpose();
    Theta_hat = X*XtXinvXt*W.asDiagonal()*Y;
    //
    for(int j = blocks(std::abs(b - 1), 0); j < blocks(std::abs(b - 1), 1) + b; j = j + 1){
      //
      Theta_tild[j] = Theta_hat[j] - (W[j]*hatvals[j]*(Y[j] - Theta_hat[j]))/(1 - W[j]*hatvals[j]);
      Tau_tild[j] = (hatvals[j] / (1 - W[j]*hatvals[j])) + pr_psi;
    }
    //
    Psi_tild[b] = pr_psi;
  }
  // Un-permute rows of Theta and Tau
  Theta_tild = perm.transpose()*Theta_tild;
  Tau_tild = perm.transpose()*Tau_tild;
  //
  Rcpp::NumericVector rTau_tild(Tau_tild.size());
  Rcpp::NumericVector rTheta_tild(Theta_tild.size());
  Rcpp::NumericVector rPsi_tild(Psi_tild.size());
  //
  for(int i = 0; i < Tau_tild.size(); i = i + 1){
      rTau_tild[i] = Tau_tild[i];
  }
  //
  for(int i = 0; i < Theta_tild.size(); i = i + 1){
      rTheta_tild[i] = Theta_tild[i];
  }
  //
  for(int i = 0; i < Psi_tild.size(); i = i + 1){
      rPsi_tild[i] = Psi_tild[i];
  }
  //
  return Rcpp::List::create(Rcpp::Named("Tau_tild") = rTau_tild,
                            Rcpp::Named("Theta_tild") = rTheta_tild,
                            Rcpp::Named("Psi_tild") = rPsi_tild);
  //
}
