#'
#'
#'
find_signal_noise_mle <- function(Y, U, V, R, sigmasq){
  # Get dimensions
  n <- nrow(U)
  m <- nrow(V)
  p <- ncol(V)
  q <- ncol(U)
  # Mean-center U, V
  if(!all(dim(U) == c(1,1))){
    U <- apply(U, 2, function(x){(x - mean(x))})
  }
  if(!all(dim(V) == c(1,1))){
    V <- apply(V, 2, function(x){(x - mean(x))})
  }
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
  Rmean <- mean(R)
  N <- prod(dim(Y))
  non_zero_svs <- svs != 0
  any_zero_svs <- any(!non_zero_svs)
  #
  opt <- function(taupsi){
    tausq <- taupsi[1]
    psisq <- taupsi[2]
    (0.5)*sum(log(((psisq)*(svs^2) + (sigmasq / Rmean) + tausq))) + # log-determinant
      (0.5)*sum((1 / ((psisq)*(svs[non_zero_svs]^2) + (sigmasq / Rmean) + tausq))*(QY^2)) + # quad form directions of X
      (0.5)*ifelse(any_zero_svs, 1, 0)*(1 / ((sigmasq / Rmean) + tausq))*(normY - normQY) + # quad form anti directions of X
      (0.95)*(normY / N)/tausq +
      (0.05)*(normY / sum(svs[non_zero_svs]^2))/psisq
  }
  #
  opt_res <- optim(c(0.1, 0.1), opt, lower = c(1e-10, 1e-10), method = "L-BFGS-B")
  opt_tausq <- opt_res[["par"]][1]
  opt_psisq <- opt_res[["par"]][2]
  return(list(tausq = opt_tausq, psisq = opt_psisq))
}
