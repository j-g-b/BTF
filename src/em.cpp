// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppDist)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <RcppDist.h>
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
Eigen::MatrixXd emRcppToEigenMat(Rcpp::NumericMatrix & RMat){
  Eigen::MatrixXd EigMat = Eigen::MatrixXd::Zero(RMat.nrow(), RMat.ncol());
  for(int i = 0; i < RMat.nrow(); i = i + 1) {
    for(int j = 0; j < RMat.ncol(); j = j + 1){
      EigMat(i, j) = RMat(i, j);
    }
  }
  return(EigMat);
}
//
Eigen::VectorXd emRcppToEigenVec(Rcpp::NumericVector & RVec){
  Eigen::VectorXd EigVec = Eigen::VectorXd::Zero(RVec.size());
  for(int i = 0; i < RVec.size(); i = i + 1) {
    EigVec[i] = RVec[i];
  }
  return(EigVec);
}
//
double emNormCDF(double x){
    return std::erfc(-x/std::sqrt(2))/2;
}
//
double ELTruncNorm(double mu, double sigmasq, double lwr){
    double alpha = (lwr - mu)/std::sqrt(sigmasq);
    double dens = (1 / std::sqrt(4*std::acos(0.0)))*std::exp(-0.5 * std::pow(alpha, 2));
    if(dens == 0){
        if(lwr < mu){
            return mu;
        } else {
            return lwr;
        }
    } else {
        double res = mu + (std::sqrt(sigmasq)*dens / (1 - emNormCDF(alpha)));
        if(isnan(res)){
            if(lwr < mu){
                return mu;
            } else {
                return lwr;
            }
        }
        if(isinf(res)){
            if(lwr < mu){
                return mu;
            } else {
                return lwr;
            }
        }
        return res;
    }
}
//
void AltMinStep(std::vector<Eigen::MatrixXd> & Tensor, std::vector<Eigen::MatrixXd> & R,
                                     Eigen::MatrixXd & U, Eigen::MatrixXd & V,
                                     Rcpp::NumericVector & MatrixType, Rcpp::List & TensorList){
    //
    Eigen::MatrixXd CU = (Tensor[0]*V)*(V.transpose()*Tensor[0].transpose());
    //
    for(int k = 1; k < MatrixType.size(); k = k + 1){
        CU = CU + (Tensor[k]*V)*(V.transpose()*Tensor[k].transpose());
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(CU);
    U = es.eigenvectors().rightCols(U.cols());
    //
    Eigen::MatrixXd CV = (Tensor[0].transpose()*U)*(U.transpose()*Tensor[0]);
    //
    for(int k = 1; k < MatrixType.size(); k = k + 1){
        CV = CV + (Tensor[k].transpose()*U)*(U.transpose()*Tensor[k]);
    }
    es.compute(CV);
    V = es.eigenvectors().rightCols(V.cols());
}
//
void MStep(std::vector<Eigen::MatrixXd> & Tensor, std::vector<Eigen::MatrixXd> & R,
                                     Eigen::MatrixXd & U, Eigen::MatrixXd & V,
                                     Rcpp::NumericVector & MatrixType, Rcpp::List & TensorList){
    //
    Eigen::MatrixXd Uprev;
    Eigen::MatrixXd Vprev;
    //
    double UMSE = 1;
    double VMSE = 1;
    //
    while(UMSE + VMSE > 1e-5){
        //
        Uprev = U;
        Vprev = V;
        //
        AltMinStep(Tensor, R, U, V, MatrixType, TensorList);
        //
        UMSE = ((U - Uprev).array().pow(2).matrix().sum()) / (U.cols()*U.rows());
        VMSE = ((V - Vprev).array().pow(2).matrix().sum()) / (V.cols()*V.rows());
        cout << UMSE + VMSE << std::endl;
    }
    //
    for(int k = 0; k < MatrixType.size(); k = k + 1){
        R[k] = U.transpose()*Tensor[k]*V;
    }
}
//
void EStep(std::vector<Eigen::MatrixXd> & Tensor, std::vector<Eigen::MatrixXd> & R,
                                     Eigen::MatrixXd & U, Eigen::MatrixXd & V,
                                     Rcpp::NumericVector & MatrixType, Rcpp::List & TensorList){
  //
  Eigen::MatrixXd PredMat;
  Rcpp::NumericMatrix TempMat;
  double SS;
  for(int k = 0; k < MatrixType.size(); k = k + 1){
    if(MatrixType[k] == 0){
      PredMat = ((U*R[k]*V.transpose()).array()).matrix();
      TempMat = Rcpp::as<NumericMatrix>(TensorList["matrix" + std::to_string(k)]);
      for(int i = 0; i < U.rows(); i = i + 1){
        for(int j = 0; j < V.rows(); j = j + 1){
          if(!isnan(TempMat(i, j))){
            Tensor[k](i, j) = TempMat(i, j);
          } else {
            Tensor[k](i, j) = PredMat(i, j);
          }
        }
      }
    }
    if(MatrixType[k] == 1){
      PredMat = ((U*R[k]*V.transpose()).array()).matrix();
      TempMat = Rcpp::as<NumericMatrix>(TensorList["matrix" + std::to_string(k)]);
      for(int i = 0; i < U.rows(); i = i + 1){
        for(int j = 0; j < V.rows(); j = j + 1){
          if(!isnan(TempMat(i, j))){
            if(TempMat(i, j) == 1){
              Tensor[k](i, j) = ELTruncNorm(PredMat(i, j), 1, 0);
            } else {
              Tensor[k](i, j) = -ELTruncNorm(-PredMat(i, j), 1, 0);
            }
          } else {
            Tensor[k](i, j) = PredMat(i ,j);
          }
        }
      }
    }
    if(MatrixType[k] == 2){
      PredMat = ((U*R[k]*V.transpose()).array()).matrix();
      TempMat = Rcpp::as<NumericMatrix>(TensorList["matrix" + std::to_string(k)]);
      SS = (Tensor[k] - PredMat).array().pow(2).matrix().sum() / (U.rows()*V.rows());
      for(int i = 0; i < U.rows(); i = i + 1){
        for(int j = 0; j < V.rows(); j = j + 1){
          if(!isnan(TempMat(i, j))){
            if(TempMat(i, j) <= 0){
              Tensor[k](i, j) = -ELTruncNorm(-PredMat(i, j), SS, 0);
            } else {
              Tensor[k](i, j) = TempMat(i, j);
            }
          } else {
            Tensor[k](i, j) = PredMat(i, j);
          }
        }
      }
    }
  }
}
//
void Impute(std::vector<Eigen::MatrixXd> & Tensor, std::vector<Eigen::MatrixXd> & R,
                                     Eigen::MatrixXd & U, Eigen::MatrixXd & V,
                                     Rcpp::NumericVector & MatrixType, Rcpp::List & TensorList){
  //
  Eigen::MatrixXd PredMat;
  for(int k = 0; k < MatrixType.size(); k = k + 1){
    PredMat = ((U*R[k]*V.transpose()).array()).matrix();
    for(int i = 0; i < U.rows(); i = i + 1){
        for(int j = 0; j < V.rows(); j = j + 1){
            Tensor[k](i, j) = PredMat(i, j);
        }
    }
  }
}
//
double Lik(std::vector<Eigen::MatrixXd> & Tensor, std::vector<Eigen::MatrixXd> & R,
                                     Eigen::MatrixXd & U, Eigen::MatrixXd & V,
                                     Rcpp::NumericVector & MatrixType, Rcpp::List & TensorList){
  //
  Eigen::MatrixXd PredMat;
  Rcpp::NumericMatrix TempMat;
  double SS;
  double LL = 0;
  for(int k = 0; k < MatrixType.size(); k = k + 1){
    if(MatrixType[k] == 0){
      PredMat = ((U*R[k]*V.transpose()).array()).matrix();
      TempMat = Rcpp::as<NumericMatrix>(TensorList["matrix" + std::to_string(k)]);
      for(int i = 0; i < U.rows(); i = i + 1){
        for(int j = 0; j < V.rows(); j = j + 1){
          if(!isnan(TempMat(i, j))){
            Tensor[k](i, j) = TempMat(i, j);
          } else {
            Tensor[k](i, j) = PredMat(i, j);
          }
        }
      }
    }
    if(MatrixType[k] == 1){
      PredMat = ((U*R[k]*V.transpose()).array()).matrix();
      TempMat = Rcpp::as<NumericMatrix>(TensorList["matrix" + std::to_string(k)]);
      for(int i = 0; i < U.rows(); i = i + 1){
        for(int j = 0; j < V.rows(); j = j + 1){
          if(!isnan(TempMat(i, j))){
            if(TempMat(i, j) == 1){
              Tensor[k](i, j) = ELTruncNorm(PredMat(i, j), 1, 0);
            } else {
              Tensor[k](i, j) = -ELTruncNorm(-PredMat(i, j), 1, 0);
            }
          } else {
            Tensor[k](i, j) = PredMat(i ,j);
          }
        }
      }
    }
    if(MatrixType[k] == 2){
      PredMat = ((U*R[k]*V.transpose()).array()).matrix();
      TempMat = Rcpp::as<NumericMatrix>(TensorList["matrix" + std::to_string(k)]);
      SS = (Tensor[k] - PredMat).array().pow(2).matrix().sum() / (U.rows()*V.rows());
      for(int i = 0; i < U.rows(); i = i + 1){
        for(int j = 0; j < V.rows(); j = j + 1){
          if(!isnan(TempMat(i, j))){
            if(TempMat(i, j) <= 0){
              Tensor[k](i, j) = -ELTruncNorm(-PredMat(i, j), SS, 0);
            } else {
              Tensor[k](i, j) = TempMat(i, j);
            }
          } else {
            Tensor[k](i, j) = PredMat(i, j);
          }
        }
      }
    }
  }
}
//
// [[Rcpp::export]]
Rcpp::List EM(Rcpp::List & TensorList, Rcpp::NumericVector & MatrixType,
               Rcpp::NumericMatrix & U0, Rcpp::NumericMatrix & V0,
               Rcpp::List & R0, Rcpp::NumericVector & SS0, Rcpp::NumericVector & Mu0,
               int S)
{
  //
  const int N = U0.nrow();
  const int M = V0.nrow();
  const int K = SS0.size();
  const int D1 = U0.ncol();
  const int D2 = V0.ncol();
  //
  // Define matrix containers for vector representations
  Eigen::MatrixXd U(N, D1);
  Eigen::MatrixXd V(M, D2);
  U = emRcppToEigenMat(U0);
  V = emRcppToEigenMat(V0);
  //
  std::vector<Eigen::MatrixXd> R(K);
  for(int k = 0; k < K; k = k + 1){
    Rcpp::NumericMatrix TempMat = Rcpp::as<NumericMatrix>(R0["matrix" + std::to_string(k)]);
    R[k] = emRcppToEigenMat(TempMat);
  }
  //
  // Extract matrices from the TensorList
  std::vector<Eigen::MatrixXd> Tensor(K);
  for(int k = 0; k < K; k = k + 1) {
      //
      Rcpp::NumericMatrix TempMat = Rcpp::as<NumericMatrix>(TensorList["matrix" + std::to_string(k)]);
      Tensor[k] = emRcppToEigenMat(TempMat);
  }
  //
  Impute(Tensor, R, U, V, MatrixType, TensorList);
  cout << "Imputing" << std::endl;
  // Run EM algorithm
  for(int s = 0; s < S; s = s + 1){
    //
    EStep(Tensor, R, U, V, MatrixType, TensorList);
    cout << "E step" << std::endl;
    //
    MStep(Tensor, R, U, V, MatrixType, TensorList);
    cout << "M step" << std::endl;
    //
    cout << s << endl;
  }
  //
  Rcpp::NumericMatrix rU(U.rows(), U.cols());
  Rcpp::NumericMatrix rV(V.rows(), V.cols());
  std::vector<Rcpp::NumericMatrix> rR(K);
  //
  for(int i = 0; i < U.rows(); i = i + 1){
    for(int j = 0; j < U.cols(); j = j + 1){
      rU(i, j) = U(i, j);
    }
  }
  //
  for(int i = 0; i < V.rows(); i = i + 1){
    for(int j = 0; j < V.cols(); j = j + 1){
      rV(i, j) = V(i, j);
    }
  }
  //
  for(int k = 0; k < R.size(); k = k + 1){
    rR[k] = Rcpp::NumericMatrix(D1, D2);
    for(int i = 0; i < D1; i = i + 1){
      for(int j = 0; j < D2; j = j + 1){
        rR[k](i, j) = R[k](i, j);
      }
    }
  }
  //
  return Rcpp::List::create(Rcpp::Named("U") = rU,
                            Rcpp::Named("V") = rV,
                            Rcpp::Named("R") = rR);
  //
}
