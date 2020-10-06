#'
#'
#'
find_signal_noise_mle <- function(Y, U, V, R, sigmasq){
  require(dfoptim)
  # Get dimensions
  n <- nrow(U)
  m <- nrow(V)
  p <- ncol(V)
  q <- ncol(U)
  # Mean-center Y
  Y <- Y - mean(Y)
  # Compute SVD of U, V
  normY <- sum(Y^2)
  if(!all(dim(U) == c(1,1))){
    svdU <- svd(U)
  }
  if(!all(dim(V) == c(1,1))){
    svdV <- svd(V)
  }
  #
  if(!all(dim(U) == c(1,1))){
    QY <- t(svdU[["u"]])%*%Y
  } else {
    QY <- Y
  }
  #
  if(!all(dim(V) == c(1,1))){
    QY <- QY%*%svdV[["u"]]
  }
  QY <- c(QY)
  normQY <- sum(QY^2)
  #
  if(!all(dim(V) == c(1,1))){
    svs <- svdV[["d"]]
    if(m > p){
      svs <- c(svs, rep(0, m - p))
    }
  } else {
    svs <- 1
  }
  #
  if(!all(dim(U) == c(1,1))){
    usvs <- svdU[["d"]]
    if(n > q){
      usvs <- c(usvs, rep(0, n - q))
    }
    svs <- c(sapply(svs, function(x){x*usvs}))
  }
  #
  Rmean <- 1/mean(1/R)
  N <- prod(dim(Y))
  non_zero_svs <- svs != 0
  any_zero_svs <- any(!non_zero_svs)
  sigmasq <- sigmasq / Inf
  #
  opt <- function(taupsi){
    tausq <- taupsi[1]
    psisq <- taupsi[2]
    (0.5)*sum(log(((psisq)*(svs^2) + (sigmasq / Rmean) + tausq))) + # log-determinant
      (0.5)*sum((1 / ((psisq)*(svs[non_zero_svs]^2) + (sigmasq / Rmean) + tausq))*(QY^2)) + # quad form directions of X
      (0.5)*ifelse(any_zero_svs, 1, 0)*(1 / ((sigmasq / Rmean) + tausq))*(normY - normQY) +# quad form anti directions of X
      (0.95)*(normY / N)/tausq +
      (0.05)*(normY / sum(svs^2))/psisq
  }
  #
  opt_res <- dfoptim::nmkb(c(0.1, 0.001), opt, lower = c(0, 0))
  opt_tausq <- opt_res[["par"]][1]
  opt_psisq <- opt_res[["par"]][2]
  #print(normQY / normY)
  #print((sigmasq / Rmean) / (normY / N))
  #print(opt_tausq*N / normY)
  #print(opt_psisq*sum(svs^2) / normY)
  return(list(tausq = opt_tausq, psisq = opt_psisq))
}
