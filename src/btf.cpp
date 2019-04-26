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
void write_to_log(std::string & text, std::string & filename){
  std::ofstream log_file(filename, std::ios_base::out | std::ios_base::app);
  log_file << text << std::endl;
}
//
arma::mat CastArma(Eigen::MatrixXd eigen_A) {
  //
  arma::mat arma_B = arma::mat(eigen_A.data(), eigen_A.rows(), eigen_A.cols());
  //
  return arma_B;
}
//
Eigen::MatrixXd CastEigen(arma::mat arma_A) {
  //
  Eigen::MatrixXd eigen_B = Eigen::Map<Eigen::MatrixXd>(arma_A.memptr(), arma_A.n_rows, arma_A.n_cols);
  //
  return eigen_B;
}
//
Eigen::MatrixXd RcppToEigenMat(Rcpp::NumericMatrix & RMat){
  Eigen::MatrixXd EigMat = Eigen::MatrixXd::Zero(RMat.nrow(), RMat.ncol());
  for(int i = 0; i < RMat.nrow(); i = i + 1) {
    for(int j = 0; j < RMat.ncol(); j = j + 1){
      EigMat(i, j) = RMat(i, j);
    }
  }
  return(EigMat);
}
//
Eigen::VectorXd RcppToEigenVec(Rcpp::NumericVector & RVec){
  Eigen::VectorXd EigVec = Eigen::VectorXd::Zero(RVec.size());
  for(int i = 0; i < RVec.size(); i = i + 1) {
    EigVec[i] = RVec[i];
  }
  return(EigVec);
}
//
double Bernoulli(double p){
  double u = arma::randu();
  if(u < p){
    return(1.0);
  } else {
    return(0.0);
  }
}
//
Eigen::MatrixXd MVNorm(int N, Eigen::VectorXd & mu, Eigen::MatrixXd & Sigma){
  //
  arma::mat aZ = arma::randn<arma::mat>(mu.size(), N);
  Eigen::MatrixXd Z = CastEigen(aZ);
  Eigen::MatrixXd L = Sigma.llt().matrixL();
  Eigen::MatrixXd X = (L*Z).colwise() + mu;
  //
  return X.transpose();
}
//
void UpdateUV(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & U, std::vector<Eigen::MatrixXd> & Tensor, std::vector<Eigen::MatrixXd> & R, 
                         Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & V, Eigen::VectorXd & MuU, Eigen::MatrixXd & LambdaU, Eigen::VectorXd & SigmaSq,
                         int NUVECS, int NVVECS, int NDIMS, int NMODS, int axis){
  //
  Eigen::VectorXd MuStar;
  Eigen::MatrixXd LambdaStarInv;
  Eigen::MatrixXd VtV = V.transpose()*V;
  if(axis == 1){
    for(int i = 0; i < NUVECS; i = i + 1){
      MuStar = Eigen::VectorXd::Zero(NDIMS);
      LambdaStarInv = Eigen::MatrixXd::Zero(NDIMS, NDIMS);
      for(int k = 0; k < NMODS; k = k + 1){
        LambdaStarInv.noalias() += (R[k]*VtV*(R[k].transpose())) / SigmaSq[k];
        MuStar.noalias() += (R[k]*(V.transpose()*(Tensor[k].row(i).asDiagonal()))).rowwise().sum() / SigmaSq[k];
      }
      LambdaStarInv = (LambdaStarInv + LambdaU).llt().solve(Eigen::MatrixXd::Identity(NDIMS, NDIMS));
      MuStar = LambdaStarInv*(MuStar + LambdaU*MuU);
      U.row(i) = MVNorm(1, MuStar, LambdaStarInv);
    }
  } else if(axis == 2) {
    for(int i = 0; i < NUVECS; i = i + 1){
      MuStar = Eigen::VectorXd::Zero(NDIMS);
      LambdaStarInv = Eigen::MatrixXd::Zero(NDIMS, NDIMS);
      for(int k = 0; k < NMODS; k = k + 1){
        LambdaStarInv.noalias() += ((R[k].transpose())*VtV*R[k]) / SigmaSq[k];
        MuStar.noalias() += ((R[k].transpose())*(V.transpose()*Tensor[k].col(i).asDiagonal())).rowwise().sum() / SigmaSq[k];
      }
      LambdaStarInv = (LambdaStarInv + LambdaU).llt().solve(Eigen::MatrixXd::Identity(NDIMS, NDIMS));
      MuStar = LambdaStarInv*(MuStar + LambdaU*MuU);
      U.row(i) = MVNorm(1, MuStar, LambdaStarInv);
    }
  }
  //
}
//
void UpdateTheta0(Eigen::VectorXd & MuU, Eigen::MatrixXd & LambdaU, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & U, Eigen::MatrixXd & W0Inv, Eigen::VectorXd & Mu0, double & beta0, double & nu0){
  //
  Eigen::VectorXd Ubar = (U.colwise().sum() / U.rows()).transpose();
  Eigen::MatrixXd S = (U.transpose()*U) / U.rows();
  //
  Eigen::VectorXd MuStar = (beta0*Mu0 + U.rows()*Ubar) / (beta0 + U.rows());
  Eigen::MatrixXd WStar = (W0Inv + U.rows()*S + ((beta0*U.rows())/(beta0 + U.rows()))*(Mu0 - Ubar)*((Mu0 - Ubar).transpose())).llt().solve(Eigen::MatrixXd::Identity(U.cols(), U.cols()));
  //
  arma::mat aWStar = CastArma(WStar);
  arma::mat aLambda0 = arma::wishrnd(aWStar, nu0 + U.rows());
  LambdaU = CastEigen(aLambda0);
  Eigen::MatrixXd LambdaStar = LambdaU.llt().solve(Eigen::MatrixXd::Identity(U.cols(), U.cols())) / (beta0 + U.rows());
  MuU = MVNorm(1, MuStar, LambdaStar).transpose();
}
//
void UpdateOmega0(Eigen::VectorXd & XiR, Eigen::MatrixXd & PsiR, std::vector<Eigen::MatrixXd> & R, Eigen::MatrixXd & Psi0Inv, Eigen::VectorXd & Xi0, double & u0, double & v0){
  //
  Eigen::VectorXd Rtemp = Eigen::VectorXd::Zero(XiR.size());
  Eigen::VectorXd Rbar = Eigen::VectorXd::Zero(XiR.size());
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(XiR.size(), XiR.size());
  for(int k = 0; k < R.size(); k = k + 1){
    Rtemp = Eigen::Map<Eigen::VectorXd>(R[k].data(), R[k].rows()*R[k].rows());
    Rbar += Rtemp;
    S += Rtemp*Rtemp.transpose();
  }
  //
  Rbar = Rbar / R.size();
  //
  Eigen::VectorXd XiStar = (u0*Xi0 + R.size()*Rbar) / (u0 + R.size());
  Eigen::MatrixXd WStar = (Psi0Inv + R.size()*S + ((u0*R.size())/(u0 + R.size()))*(Xi0 - Rbar)*((Xi0 - Rbar).transpose())).llt().solve(Eigen::MatrixXd::Identity(XiR.size(), XiR.size()));
  //
  arma::mat aWStar = CastArma(WStar);
  arma::mat aPsi0 = arma::wishrnd(aWStar, v0 + R.size());
  PsiR = CastEigen(aPsi0);
  Eigen::MatrixXd PsiStar = PsiR.llt().solve(Eigen::MatrixXd::Identity(XiR.size(), XiR.size())) / (u0 + R.size());
  XiR = MVNorm(1, XiStar, PsiStar).transpose();
}
//
void UpdateSigmaSq(Eigen::VectorXd & SigmaSq, std::vector<Eigen::MatrixXd> & Tensor, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & U, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & V, std::vector<Eigen::MatrixXd> & R, Rcpp::NumericVector & MatrixType, double n0, double eta0){
  //
  double SS;
  double post_a;
  double post_b;
  double ss;
  double sample_size;
  //
  for(int k = 0; k < SigmaSq.size(); k = k + 1){
    if(MatrixType[k] == 0){
      SS = 0;
      sample_size = 0;
      for(int i = 0; i < Tensor[k].rows(); i = i + 1){
        for(int j = 0; j < Tensor[k].cols(); j = j + 1){
          if(!isnan(Tensor[k](i, j))){
            ss = Tensor[k](i, j) - U.row(i)*R[k]*(V.row(j).transpose());
            SS += pow(ss, 2);
            sample_size += 1;
          }
        }
      }
      post_a = (n0 + sample_size) / 2.0;
      post_b = (n0*eta0 + SS) / 2.0;
      SigmaSq[k] = 1.0 / (arma::randg(arma::distr_param(post_a, 1.0 / post_b)));
    } else {
      SigmaSq[k] = 1.0;
    }
  }
}
//
void UpdateR(std::vector<Eigen::MatrixXd> & R, std::vector<Eigen::MatrixXd> & Tensor, 
                                     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & U, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & V, Eigen::VectorXd & SigmaSq, Eigen::VectorXd & Beta0, Eigen::MatrixXd & Sigma0){
  //
  Eigen::MatrixXd SigmaStarInv;
  Eigen::VectorXd BetaStar;
  Eigen::MatrixXd BetaStarMat;
  Eigen::MatrixXd TempMat;
  Eigen::VectorXd TempVec;
  Eigen::MatrixXd UtU;
  Eigen::MatrixXd VtV;
  //
  for(int k = 0; k < R.size(); k = k + 1){
    //
    UtU = Eigen::MatrixXd::Zero(U.cols(), U.cols());
    VtV = Eigen::MatrixXd::Zero(V.cols(), V.cols());
    SigmaStarInv = Eigen::MatrixXd::Zero(U.cols()*U.cols(), U.cols()*U.cols());
    //
    BetaStar = Eigen::VectorXd::Zero(U.cols()*U.cols());
    BetaStarMat = U.transpose()*Tensor[k]*V;
    //
    UtU = U.transpose()*U;
    VtV = V.transpose()*V;
    for(int d1 = 0; d1 < V.cols(); d1 = d1 + 1){
      for(int d2 = d1; d2 < V.cols(); d2 = d2 + 1){
        SigmaStarInv.block(d1*V.cols(), d2*V.cols(), V.cols(), V.cols()).noalias() += VtV(d1, d2)*UtU;
      }
    }
    //
    for(int d1 = 0; d1 < V.cols()*V.cols(); d1 = d1 + 1){
      for(int d2 = d1; d2 < V.cols()*V.cols(); d2 = d2 + 1){
        SigmaStarInv(d2, d1) = SigmaStarInv(d1, d2);
      }
    }
    //
    SigmaStarInv = ((SigmaStarInv / SigmaSq[k]) + Sigma0).llt().solve(Eigen::MatrixXd::Identity(U.cols()*U.cols(), U.cols()*U.cols()));
    //
    BetaStar = Eigen::Map<Eigen::VectorXd>(BetaStarMat.data(), U.cols()*U.cols());
    BetaStar = SigmaStarInv*((BetaStar / SigmaSq[k]) + Sigma0*Beta0);
    //
    BetaStar = MVNorm(1, BetaStar, SigmaStarInv).transpose();
    //
    R[k] = Eigen::Map<Eigen::MatrixXd>(BetaStar.data(), V.cols(), V.cols());
  }
}
//
double EvalMSE(std::vector<Eigen::MatrixXd> & Tensor, std::vector<Eigen::MatrixXd> & R, Rcpp::NumericVector & MatrixType,
                                     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & U, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & V){
  double MSE = 0;
  int counter = 0;
  for(int k = 0; k < Tensor.size(); k = k + 1) {
    if(MatrixType[k] != 1){
      for(int i = 0; i < Tensor[k].rows(); i = i + 1) {
        for(int j = 0; j < Tensor[k].cols(); j = j + 1){
          if(!isnan(Tensor[k](i, j))){
            MSE += pow(Tensor[k](i, j) - U.row(i)*R[k]*(V.row(j).transpose()), 2);
            counter += 1;
          }
        }
      }
    }
  }
  return(MSE / counter);
}
//
// Inverse cumulative distribution function (aka the probit function)
double iNormCDF(double quantile) {
  // This is the Beasley-Springer-Moro algorithm which can 
  // be found in Glasserman [2004].
  double a[4] = {   2.50662823884,
                         -18.61500062529,
                          41.39119773534,
                         -25.44106049637};

  double b[4] = {  -8.47351093090,
                          23.08336743743,
                         -21.06224101826,
                           3.13082909833};

  double c[9] = {0.3374754822726147,
                        0.9761690190917186,
                        0.1607979714918209,
                        0.0276438810333863,
                        0.0038405729373609,
                        0.0003951896511919,
                        0.0000321767881768,
                        0.0000002888167364,
                        0.0000003960315187};

  if (quantile >= 0.5 && quantile <= 0.92) {
    double num = 0.0;
    double denom = 1.0;

    for (int i=0; i<4; i++) {
      num += a[i] * pow((quantile - 0.5), 2*i + 1);
      denom += b[i] * pow((quantile - 0.5), 2*i);
    }
    return num/denom;

  } else if (quantile > 0.92 && quantile < 1) {
    double num = 0.0;

    for (int i=0; i<9; i++) {
      num += c[i] * pow((log(-log(1-quantile))), i);
    }
    return num;

  } else {
    return -1.0*iNormCDF(1-quantile);
  }
}
//
double LTruncNorm(double mu, double sigmasq, double lwr){
  double u = arma::randu();
  double z = std::sqrt(sigmasq)*iNormCDF(u) + mu;
  while(z < lwr){
    u = arma::randu();
    z = std::sqrt(sigmasq)*iNormCDF(u) + mu;
  }
  return(z);
}
//
double UTruncNorm(double mu, double sigmasq, double upr){
  double u = arma::randu();
  double z = std::sqrt(sigmasq)*iNormCDF(u) + mu;
  while(z > upr){
    u = arma::randu();
    z = std::sqrt(sigmasq)*iNormCDF(u) + mu;
  }
  return(z);
}
//
double NormCDF(double x) // Phi(-∞, x) aka N(x)
{
    return std::erfc(-x/std::sqrt(2))/2;
}
//
void UpdateZ(std::vector<Eigen::MatrixXd> & Tensor, std::vector<Eigen::MatrixXd> & R, 
                                     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & U, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & V,
                                     Rcpp::NumericVector & MatrixType, Rcpp::List & TensorList){
  //
  Eigen::MatrixXd PredMat;
  Rcpp::NumericMatrix TempMat;
  for(int k = 0; k < MatrixType.size(); k = k + 1){
    if(MatrixType[k] == 1){
      PredMat = U*R[k]*V.transpose();
      TempMat = Rcpp::as<NumericMatrix>(TensorList["matrix" + std::to_string(k)]);
      for(int i = 0; i < U.rows(); i = i + 1){
        for(int j = 0; j < V.rows(); j = j + 1){
          if(!isnan(TempMat(i, j))){
            if(TempMat(i, j) == 1){
              Tensor[k](i, j) = LTruncNorm(PredMat(i, j), 1, 0);
            } else {
              Tensor[k](i, j) = UTruncNorm(PredMat(i, j), 1, 0);
            }
          }
        }
      }
    }
  }
}
//
void UpdateMissing(std::vector<Eigen::MatrixXd> & Tensor, std::vector<Eigen::MatrixXd> & R, 
                                     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & U, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & V,
                                     Eigen::VectorXd & SigmaSq, Rcpp::NumericVector & MatrixType, Rcpp::List & TensorList){
  //
  Eigen::MatrixXd PredMat;
  Rcpp::NumericMatrix TempMat;
  for(int k = 0; k < MatrixType.size(); k = k + 1){
    if(MatrixType[k] == 0){
      PredMat = U*R[k]*V.transpose();
      TempMat = Rcpp::as<NumericMatrix>(TensorList["matrix" + std::to_string(k)]);
      double sigma = std::sqrt(SigmaSq[k]);
      for(int i = 0; i < U.rows(); i = i + 1){
        for(int j = 0; j < V.rows(); j = j + 1){
          if(isnan(TempMat(i, j))){
            Tensor[k](i, j) = PredMat(i, j) + sigma*arma::randn();
          }
        }
      }
    } else {
      PredMat = U*R[k]*V.transpose();
      TempMat = Rcpp::as<NumericMatrix>(TensorList["matrix" + std::to_string(k)]);
      for(int i = 0; i < U.rows(); i = i + 1){
        for(int j = 0; j < V.rows(); j = j + 1){
          if(isnan(TempMat(i, j))){
            Tensor[k](i, j) = Bernoulli(NormCDF(PredMat(i, j)));
          }
        }
      }
    }
  }
}
//
void WriteMatCSV(std::string name, Eigen::MatrixXd matrix)
{
  std::ofstream file(name.c_str());

  for(int  i = 0; i < matrix.rows(); i++){
      for(int j = 0; j < matrix.cols(); j++){
         std::string str = std::to_string(matrix(i,j));
         if(j+1 == matrix.cols()){
             file<<str;
         }else{
             file<<str<<',';
         }
      }
      file<<'\n';
  }
  file.close();
}
//
void WriteVecCSV(std::string name, Eigen::VectorXd vec)
{
  std::ofstream file(name.c_str());

  for(int  i = 0; i < vec.size(); i++){
    std::string str = std::to_string(vec[i]);
    file<<str<<'\n';
  }
  file.close();
}
//
void SaveState(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & U, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & V, 
               std::vector<Eigen::MatrixXd> & R, Eigen::VectorXd & SigmaSq, 
               Eigen::MatrixXd & LambdaU, Eigen::MatrixXd & LambdaV, 
               Eigen::VectorXd & MuU, Eigen::VectorXd & MuV, 
               Eigen::VectorXd & XiR, Eigen::MatrixXd & PsiR, 
               int S, std::string dir_name){
  //
  std::string fstring;
  //
  fstring = dir_name + "/U_" + std::to_string(S) + ".csv";
  WriteMatCSV(fstring, U);
  fstring = dir_name + "/V_" + std::to_string(S) + ".csv";
  WriteMatCSV(fstring, V);
  for(int k = 0; k < R.size(); k = k + 1){
    fstring = dir_name + "/R" + std::to_string(k) + "_" + std::to_string(S) + ".csv";
    WriteMatCSV(fstring, R[k]);
  }
  fstring = dir_name + "/SS_" + std::to_string(S) + ".csv";
  WriteVecCSV(fstring, SigmaSq);
  fstring = dir_name + "/LambdaU_" + std::to_string(S) + ".csv";
  WriteMatCSV(fstring, LambdaU);
  fstring = dir_name + "/LambdaV_" + std::to_string(S) + ".csv";
  WriteMatCSV(fstring, LambdaV);
  fstring = dir_name + "/MuU_" + std::to_string(S) + ".csv";
  WriteVecCSV(fstring, MuU);
  fstring = dir_name + "/MuV_" + std::to_string(S) + ".csv";
  WriteVecCSV(fstring, MuV);
  fstring = dir_name + "/XiR_" + std::to_string(S) + ".csv";
  WriteVecCSV(fstring, XiR);
    fstring = dir_name + "/PsiR_" + std::to_string(S) + ".csv";
  WriteMatCSV(fstring, PsiR);
  //
}
//
// [[Rcpp::export]]
Rcpp::List BTF(Rcpp::List & TensorList, Rcpp::NumericVector & MatrixType, 
               Rcpp::NumericMatrix & U0, Rcpp::NumericMatrix & V0, 
               Rcpp::List & R0, Rcpp::NumericVector & SS0, 
               Rcpp::NumericMatrix & LambdaU0, Rcpp::NumericMatrix & LambdaV0, 
               Rcpp::NumericVector & MuU0, Rcpp::NumericVector & MuV0, 
               Rcpp::NumericVector & Xi0, Rcpp::NumericMatrix & Psi0, 
               int S, SEXP SaveDir, int Thin, int Burn, int Start)
{
  //
  std::string save_dir = Rcpp::as<std::string>(SaveDir);
  //
  const int N = U0.nrow();
  const int M = V0.nrow();
  const int K = SS0.size();
  const int D = U0.ncol();
  //
  // Define matrix containers for vector representations
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> U(N, D);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> V(M, D);
  U = RcppToEigenMat(U0);
  V = RcppToEigenMat(V0);
  //
  std::vector<Eigen::MatrixXd> R(K);
  for(int k = 0; k < K; k = k + 1){
    Rcpp::NumericMatrix TempMat = Rcpp::as<NumericMatrix>(R0["matrix" + std::to_string(k)]);
    R[k] = RcppToEigenMat(TempMat);
  }
  //
  Eigen::VectorXd SigmaSq = RcppToEigenVec(SS0);
  // Extract matrices from the TensorList
  std::vector<Eigen::MatrixXd> Tensor(K);
  for(int k = 0; k < K; k = k + 1) {
      //
      Rcpp::NumericMatrix TempMat = Rcpp::as<NumericMatrix>(TensorList["matrix" + std::to_string(k)]);
      Tensor[k] = RcppToEigenMat(TempMat);
  }
  //
  Eigen::VectorXd MuU = RcppToEigenVec(MuU0);
  Eigen::VectorXd MuV = RcppToEigenVec(MuV0);
  Eigen::MatrixXd LambdaU = RcppToEigenMat(LambdaU0);
  Eigen::MatrixXd LambdaV = RcppToEigenMat(LambdaV0);
  Eigen::VectorXd XiR = RcppToEigenVec(Xi0);
  Eigen::MatrixXd PsiR = RcppToEigenMat(Psi0);
  //
  // Set up containers for parameters to return
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Umean(N, D);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Vmean(M, D);
  Eigen::VectorXd SSmean = Eigen::VectorXd::Zero(K);
  std::vector<Eigen::MatrixXd> Rmean(K);
  for(int k = 0; k < K; k = k + 1){
    Rmean[k] = Eigen::MatrixXd::Zero(D, D);
  }
  //
  // Set initial parameters for Gibbs
  Eigen::VectorXd Mu0 = Eigen::VectorXd::Zero(D);
  Eigen::MatrixXd W0Inv = Eigen::MatrixXd::Identity(D, D);
  Eigen::VectorXd XiR0 = Eigen::VectorXd::Zero(D*D);
  Eigen::MatrixXd PsiR0 = Eigen::MatrixXd::Identity(D*D, D*D);
  double beta0 = 1;
  double nu0 = D;
  double n0 = 5;
  double eta0 = 5;
  double u0 = 1;
  double v0 = D*D;
  //
  // Run Gibbs sampler
  for(int s = 0; s < S; s = s + 1){
    //
    UpdateMissing(Tensor, R, U, V, SigmaSq, MatrixType, TensorList);
    cout << "Update Missing" << endl;
    //
    UpdateZ(Tensor, R, U, V, MatrixType, TensorList);
    cout << "Update Z" << endl;
    //
    UpdateTheta0(MuU, LambdaU, U, W0Inv, Mu0, beta0, nu0);
    UpdateTheta0(MuV, LambdaV, V, W0Inv, Mu0, beta0, nu0);
    cout << "Update Theta" << endl;
    UpdateOmega0(XiR, PsiR, R, PsiR0, XiR0, u0, v0);
    cout << "Update Omega" << endl;
    UpdateUV(U, Tensor, R, V, MuU, LambdaU, SigmaSq, N, M, D, K, 1);
    cout << "Update U" << endl;
    UpdateUV(V, Tensor, R, U, MuV, LambdaV, SigmaSq, M, N, D, K, 2);
    cout << "Update V" << endl;
    UpdateR(R, Tensor, U, V, SigmaSq, XiR, PsiR);
    cout << "Update R" << endl;
    UpdateSigmaSq(SigmaSq, Tensor, U, V, R, MatrixType, n0, eta0);
    //
    if(s > Burn - 1){
      //
      // cout << EvalMSE(Tensor, Rmean, MatrixType, Umean, Vmean) << endl;
      if(s % Thin == 0){
        if(save_dir != "0"){
          SaveState(U, V, R, SigmaSq, LambdaU, LambdaV, MuU, MuV, XiR, PsiR, s + Start, save_dir);
        }
      }
    }
  }
  Umean = U;
  Vmean = V;
  SSmean = SigmaSq;
  for(int k = 0; k < K; k = k + 1){
    Rmean[k] = R[k];
  }
  //
  Rcpp::NumericMatrix rU(Umean.rows(), Umean.cols());
  Rcpp::NumericMatrix rV(Vmean.rows(), Vmean.cols());
  Rcpp::NumericVector rSS(SigmaSq.size());
  std::vector<Rcpp::NumericMatrix> rR(K);
  //
  for(int i = 0; i < Umean.rows(); i = i + 1){
    for(int j = 0; j < Umean.cols(); j = j + 1){
      rU(i, j) = Umean(i, j);
    }
  }
  //
  for(int i = 0; i < Vmean.rows(); i = i + 1){
    for(int j = 0; j < Vmean.cols(); j = j + 1){
      rV(i, j) = Vmean(i, j);
    }
  }
  //
  for(int k = 0; k < SigmaSq.size(); k = k + 1){
    rSS[k] = SigmaSq[k];
  }
  //
  for(int k = 0; k < R.size(); k = k + 1){
    rR[k] = Rcpp::NumericMatrix(D, D);
    for(int i = 0; i < D; i = i + 1){
      for(int j = 0; j < D; j = j + 1){
        rR[k](i, j) = Rmean[k](i, j);
      }
    }
  }
  //
  return Rcpp::List::create(Rcpp::Named("U") = rU,
                            Rcpp::Named("V") = rV,
                            Rcpp::Named("SS") = rSS,
                            Rcpp::Named("R") = rR);
  //
}
