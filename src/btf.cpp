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
// Input mu matrix has different means as columns, dimension of mvnorm as rows
Eigen::MatrixXd MatMVNorm(Eigen::MatrixXd & mu, Eigen::MatrixXd & Sigma){
  //
  arma::mat aZ = arma::randn<arma::mat>(mu.rows(), mu.cols());
  Eigen::MatrixXd Z = CastEigen(aZ);
  Eigen::MatrixXd L = Sigma.llt().matrixL();
  Eigen::MatrixXd X = L*Z + mu;
  //
  return X.transpose();
}
//
void UpdateUV(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & U, std::vector<Eigen::MatrixXd> & Tensor, std::vector<Eigen::MatrixXd> & R,
                         Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & V, Eigen::VectorXd & SigmaSq, Eigen::VectorXd & Mu,
                         int NUVECS, int NVVECS, int NDIMS, int NMODS, int axis){
  //
  Eigen::MatrixXd MuStar;
  Eigen::MatrixXd LambdaStarInv;
  Eigen::MatrixXd VtV = V.transpose()*V;
  if(axis == 1){
    MuStar = Eigen::MatrixXd::Zero(NDIMS, NUVECS);
    LambdaStarInv = Eigen::MatrixXd::Zero(NDIMS, NDIMS);
    for(int k = 0; k < NMODS; k = k + 1){
      LambdaStarInv.noalias() += (R[k]*VtV*(R[k].transpose())) / SigmaSq[k];
      MuStar.noalias() += (R[k]*(V.transpose())*(Tensor[k].transpose().array() - Mu[k]).matrix()) / SigmaSq[k];
    }
    LambdaStarInv = (LambdaStarInv + Eigen::MatrixXd::Identity(NDIMS, NDIMS)).llt().solve(Eigen::MatrixXd::Identity(NDIMS, NDIMS));
    MuStar = LambdaStarInv*MuStar;
    U = MatMVNorm(MuStar, LambdaStarInv);
  } else if(axis == 2) {
    MuStar = Eigen::MatrixXd::Zero(NDIMS, NUVECS);
    LambdaStarInv = Eigen::MatrixXd::Zero(NDIMS, NDIMS);
    for(int k = 0; k < NMODS; k = k + 1){
      LambdaStarInv.noalias() += (R[k].transpose()*VtV*(R[k])) / SigmaSq[k];
      MuStar.noalias() += (R[k].transpose()*(V.transpose())*(Tensor[k].array() - Mu[k]).matrix()) / SigmaSq[k];
    }
    LambdaStarInv = (LambdaStarInv + Eigen::MatrixXd::Identity(NDIMS, NDIMS)).llt().solve(Eigen::MatrixXd::Identity(NDIMS, NDIMS));
    MuStar = LambdaStarInv*MuStar;
    U = MatMVNorm(MuStar, LambdaStarInv);
  }
  //
}
//
void UpdateSigmaSq(Eigen::VectorXd & SigmaSq, std::vector<Eigen::MatrixXd> & Tensor, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & U, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & V, std::vector<Eigen::MatrixXd> & R, Eigen::VectorXd & Mu, Rcpp::NumericVector & MatrixType){
  //
  double SS;
  double post_a;
  double post_b;
  double sample_size;
  //
  for(int k = 0; k < SigmaSq.size(); k = k + 1){
    if(MatrixType[k] == 1) {
      SigmaSq[k] = 1.0;
    } else {
      SS = ((Tensor[k] - U*R[k]*V.transpose()).array() - Mu[k]).square().matrix().sum();
      sample_size = U.rows()*V.rows();
      post_a = (1 + sample_size) / 2.0;
      post_b = (1 + SS) / 2.0;
      SigmaSq[k] = 1.0 / (arma::randg(arma::distr_param(post_a, 1.0 / post_b)));
    }
  }
}
//
void UpdateMu(Eigen::VectorXd & Mu, std::vector<Eigen::MatrixXd> & Tensor, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & U, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & V, std::vector<Eigen::MatrixXd> & R, Eigen::VectorXd & SigmaSq, Rcpp::NumericVector & MatrixType){
  //
  double Gam;
  double Phi;
  //
  for(int k = 0; k < SigmaSq.size(); k = k + 1){
    Phi = ((U.rows()*V.rows()) / SigmaSq[k]) + 1;
    Gam = (Tensor[k] - U*R[k]*V.transpose()).sum() / (Phi*SigmaSq[k]);
    Mu[k] = Gam + std::sqrt(1 / Phi)*arma::randn();
  }
}
//
void UpdateR(std::vector<Eigen::MatrixXd> & R, std::vector<Eigen::MatrixXd> & Tensor,
                                     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & U, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & V, Eigen::VectorXd & SigmaSq, Eigen::VectorXd & Mu){
  //
  Eigen::MatrixXd BetaStarMat;
  Eigen::MatrixXd UtU = U.transpose()*U;
  Eigen::MatrixXd VtV = V.transpose()*V;
  Eigen::MatrixXd UtUinv = UtU.llt().solve(Eigen::MatrixXd::Identity(U.cols(), U.cols()));
  Eigen::MatrixXd VtVinv = VtV.llt().solve(Eigen::MatrixXd::Identity(V.cols(), V.cols()));
  Eigen::MatrixXd LUtUinv = UtUinv.llt().matrixL();
  Eigen::MatrixXd LVtVinv = VtVinv.llt().matrixL();
  arma::mat aZ = arma::randn<arma::mat>(U.cols(), V.cols());
  Eigen::MatrixXd Z = CastEigen(aZ);
  //
  for(int k = 0; k < R.size(); k = k + 1){
    //
    BetaStarMat = UtUinv*U.transpose()*(Tensor[k].array() - Mu[k]).matrix()*V*VtVinv;
    //
    R[k] = BetaStarMat + (std::sqrt(SigmaSq[k])*LUtUinv*Z*LVtVinv.transpose());
    //
    aZ = arma::randn<arma::mat>(U.cols(), V.cols());
    Z = CastEigen(aZ);
  }
}
//
double EvalMSE(std::vector<Eigen::MatrixXd> & Tensor, std::vector<Eigen::MatrixXd> & R, Rcpp::NumericVector & MatrixType,
                                     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & U, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & V, Eigen::VectorXd & Mu){
  double MSE = 0;
  int counter = 0;
  for(int k = 0; k < Tensor.size(); k = k + 1) {
    if(MatrixType[k] != 1){
      for(int i = 0; i < Tensor[k].rows(); i = i + 1) {
        for(int j = 0; j < Tensor[k].cols(); j = j + 1){
          if(!isnan(Tensor[k](i, j))){
            MSE += pow(Tensor[k](i, j) - U.row(i)*R[k]*(V.row(j).transpose()) - Mu[k], 2);
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
double NormCDF(double x){
    return std::erfc(-x/std::sqrt(2))/2;
}
//
double LTruncNorm(double mu, double sigmasq, double lwr){
  double z;
  z = r_truncnorm(mu, std::sqrt(sigmasq), lwr, std::numeric_limits<double>::infinity());
  if(isinf(z)){
    z = lwr;
  }
  return z;
}
//
void UpdateZ(std::vector<Eigen::MatrixXd> & Tensor, std::vector<Eigen::MatrixXd> & R,
                                     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & U, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & V,
                                     Eigen::VectorXd & Mu, Eigen::VectorXd & SigmaSq,
                                     Rcpp::NumericVector & MatrixType, Rcpp::List & TensorList){
  //
  Eigen::MatrixXd PredMat;
  Rcpp::NumericMatrix TempMat;
  for(int k = 0; k < MatrixType.size(); k = k + 1){
    if(MatrixType[k] == 1){
      PredMat = ((U*R[k]*V.transpose()).array() + Mu[k]).matrix();
      TempMat = Rcpp::as<NumericMatrix>(TensorList["matrix" + std::to_string(k)]);
      for(int i = 0; i < U.rows(); i = i + 1){
        for(int j = 0; j < V.rows(); j = j + 1){
          if(!isnan(TempMat(i, j))){
            if(TempMat(i, j) == 1){
              Tensor[k](i, j) = LTruncNorm(PredMat(i, j), 1, 0);
            } else {
              Tensor[k](i, j) = -LTruncNorm(-PredMat(i, j), 1, 0);
            }
          } else {
            if(Tensor[k](i, j) == 1){
              Tensor[k](i, j) = LTruncNorm(PredMat(i, j), 1, 0);
            } else {
              Tensor[k](i, j) = -LTruncNorm(-PredMat(i, j), 1, 0);
            }
          }
        }
      }
    }
    if(MatrixType[k] == 2){
      PredMat = ((U*R[k]*V.transpose()).array() + Mu[k]).matrix();
      TempMat = Rcpp::as<NumericMatrix>(TensorList["matrix" + std::to_string(k)]);
      for(int i = 0; i < U.rows(); i = i + 1){
        for(int j = 0; j < V.rows(); j = j + 1){
          if(!isnan(TempMat(i, j))){
            if(TempMat(i, j) <= 0){
              Tensor[k](i, j) = -LTruncNorm(-PredMat(i, j), SigmaSq[k], 0);
            }
          } else {
            if(Tensor[k](i, j) <= 0){
              Tensor[k](i, j) = -LTruncNorm(-PredMat(i, j), SigmaSq[k], 0);
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
                                     Eigen::VectorXd & SigmaSq, Eigen::VectorXd & Mu, Rcpp::NumericVector & MatrixType, Rcpp::List & TensorList){
  //
  Eigen::MatrixXd PredMat;
  Rcpp::NumericMatrix TempMat;
  for(int k = 0; k < MatrixType.size(); k = k + 1){
    if(MatrixType[k] == 0){
      PredMat = ((U*R[k]*V.transpose()).array() + Mu[k]).matrix();
      TempMat = Rcpp::as<NumericMatrix>(TensorList["matrix" + std::to_string(k)]);
      double sigma = std::sqrt(SigmaSq[k]);
      for(int i = 0; i < U.rows(); i = i + 1){
        for(int j = 0; j < V.rows(); j = j + 1){
          if(isnan(TempMat(i, j))){
            Tensor[k](i, j) = PredMat(i, j) + sigma*arma::randn();
          }
        }
      }
    }
    if(MatrixType[k] == 1){
      PredMat = ((U*R[k]*V.transpose()).array() + Mu[k]).matrix();
      TempMat = Rcpp::as<NumericMatrix>(TensorList["matrix" + std::to_string(k)]);
      for(int i = 0; i < U.rows(); i = i + 1){
        for(int j = 0; j < V.rows(); j = j + 1){
          if(isnan(TempMat(i, j))){
            Tensor[k](i, j) = Bernoulli(NormCDF(PredMat(i, j)));
          }
        }
      }
    }
    if(MatrixType[k] == 2){
      PredMat = ((U*R[k]*V.transpose()).array() + Mu[k]).matrix();
      TempMat = Rcpp::as<NumericMatrix>(TensorList["matrix" + std::to_string(k)]);
      double sigma = std::sqrt(SigmaSq[k]);
      for(int i = 0; i < U.rows(); i = i + 1){
        for(int j = 0; j < V.rows(); j = j + 1){
          if(isnan(TempMat(i, j))){
            Tensor[k](i, j) = PredMat(i, j) + sigma*arma::randn();
            if(Tensor[k](i, j) <= 0){
              Tensor[k](i, j) = 0;
            }
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

  for(int i = 0; i < matrix.rows(); i++){
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
               std::vector<Eigen::MatrixXd> & R, Eigen::VectorXd & SigmaSq, Eigen::VectorXd & Mu,
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
  fstring = dir_name + "/Mu_" + std::to_string(S) + ".csv";
  WriteVecCSV(fstring, Mu);
  //
}
//
// [[Rcpp::export]]
Rcpp::List BTF(Rcpp::List & TensorList, Rcpp::NumericVector & MatrixType,
               Rcpp::NumericMatrix & U0, Rcpp::NumericMatrix & V0,
               Rcpp::List & R0, Rcpp::NumericVector & SS0, Rcpp::NumericVector & Mu0,
               int S, SEXP SaveDir, int Thin, int Burn, int Start)
{
  //
  std::string save_dir = Rcpp::as<std::string>(SaveDir);
  //
  const int N = U0.nrow();
  const int M = V0.nrow();
  const int K = SS0.size();
  const int D1 = U0.ncol();
  const int D2 = V0.ncol();
  //
  // Define matrix containers for vector representations
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> U(N, D1);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> V(M, D2);
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
  Eigen::VectorXd Mu = RcppToEigenVec(Mu0);
  // Extract matrices from the TensorList
  std::vector<Eigen::MatrixXd> Tensor(K);
  for(int k = 0; k < K; k = k + 1) {
      //
      Rcpp::NumericMatrix TempMat = Rcpp::as<NumericMatrix>(TensorList["matrix" + std::to_string(k)]);
      Tensor[k] = RcppToEigenMat(TempMat);
  }
  //
  // Set up containers for parameters to return
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Umean(N, D1);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Vmean(M, D2);
  Eigen::VectorXd SSmean = Eigen::VectorXd::Zero(K);
  std::vector<Eigen::MatrixXd> Rmean(K);
  for(int k = 0; k < K; k = k + 1){
    Rmean[k] = Eigen::MatrixXd::Zero(D1, D2);
  }
  //
  // Run Gibbs sampler
  for(int s = 0; s < S; s = s + 1){
    //
    UpdateMissing(Tensor, R, U, V, SigmaSq, Mu, MatrixType, TensorList);
    //cout << "Update Missing" << endl;
    //
    UpdateZ(Tensor, R, U, V, Mu, SigmaSq, MatrixType, TensorList);
    //cout << "Update Z" << endl;
    //
    UpdateUV(U, Tensor, R, V, SigmaSq, Mu, N, M, D1, K, 1);
    //cout << "Update U" << endl;
    UpdateUV(V, Tensor, R, U, SigmaSq, Mu, M, N, D2, K, 2);
    //cout << "Update V" << endl;
    UpdateR(R, Tensor, U, V, SigmaSq, Mu);
    //cout << "Update R" << endl;
    UpdateSigmaSq(SigmaSq, Tensor, U, V, R, Mu, MatrixType);
    UpdateMu(Mu, Tensor, U, V, R, SigmaSq, MatrixType);
    //
    cout << s << endl;
    if(s > Burn - 1){
      //
      // cout << EvalMSE(Tensor, R, MatrixType, U, V, Mu) << endl;
      if(s % Thin == 0){
        // cout << EvalMSE(Tensor, R, MatrixType, U, V, Mu) << endl;
        cout << "Gibbs iterate: " + std::to_string(s) + "\r";
        if(save_dir != "0"){
          SaveState(U, V, R, SigmaSq, Mu, s + Start, save_dir);
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
    rR[k] = Rcpp::NumericMatrix(D1, D2);
    for(int i = 0; i < D1; i = i + 1){
      for(int j = 0; j < D2; j = j + 1){
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
