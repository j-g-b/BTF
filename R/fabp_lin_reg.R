#
#' Compute FAB p-values under linear regression linking model
#'
#' @description Comutes FAB p-values based on a linear regression linking model for matrix data given row and column covariates.
#'
#' @param Y matrix of averages of experimental readout values over R replicates
#' @param S matrix of standard errors of experimental readout values over R replicates
#' @param R number of replicates per readout value
#' @param U the row features
#' @param V the column features
#' @param pool_sampling_var logical; indicates whether sampling variance should be assumed the same across hypothesis tests
#' @param Y1 matrix of averages of experimental readout values over R1 replicates (used for contrast scores for two-sample t-tests; Y is taken to be the other sample)
#' @param S1 matrix of standard errors of experimental readout values in Y1 over R1 replicates
#' @param R1 number of replicates per readout value in Y1
#'
#' @return A data.frame of FAB p-values and the standard UMP p-values, one for each entry in Y (or each contrast score Y - Y1)
#'
#' @export fabp_lin_reg
#'
fabp_lin_reg <- function(Y, S, R, U, V, pool_sampling_var = T, Y1 = NULL, S1 = NULL, R1 = NULL){
  # Get experimental dimensions
  vecY <- c(Y)
  vecR <- c(R)
  vecS <- c(S)
  m <- nrow(V)
  n <- nrow(U)
  p <- ncol(V)
  q <- ncol(U)
  # Split data and estimate sigmasq_hat and sigmasq_tild
  hat_indices <- sample(m*n, round(m*n / 2))
  tild_indices <- setdiff(1:(m*n), hat_indices)
  if(is.null(Y1)){
    hat_nu <- sum(vecR[hat_indices] - 1)
    tild_nu <- sum(vecR[tild_indices] - 1)
    sigmasq_hat <- sum((vecR[hat_indices] - 1)*vecS[hat_indices]^2) / hat_nu
    sigmasq_tild <- sum((vecR[tild_indices] - 1)*vecS[tild_indices]^2) / tild_nu
  } else {
    hat_nu <- sum(vecR[hat_indices] - 1) + sum(c(R1)[hat_indices] - 1)
    tild_nu <- sum(vecR[tild_indices] - 1) + sum(c(R1)[tild_indices] - 1)
    sigmasq_hat <- sum((vecR[hat_indices] - 1)*(vecS[hat_indices]^2) + (c(R1)[hat_indices] - 1)*((c(S1)[hat_indices])^2)) / hat_nu
    sigmasq_tild <- sum((vecR[tild_indices] - 1)*(vecS[tild_indices]^2) + (c(R1)[tild_indices] - 1)*((c(S1)[tild_indices])^2)) / tild_nu
    vecY <- vecY - c(Y1)
    Y <- Y - Y1
    R <- 1/((1/R) + (1/R1))
  }
  # Split data to estimate signal to noise
  larger_dim <- ifelse(nrow(Y) >= ncol(Y), "row", "col")
  if(larger_dim == "row"){
    #
    rand_rows <- sample(1:nrow(Y), round(nrow(Y)/2))
    mle_taupsi <- BTF::find_signal_noise_mle(Y[rand_rows, , drop = F], U[rand_rows, ], V, R[rand_rows, , drop = F], sigmasq_tild)
    opt_tausq1 <- mle_taupsi[["tausq"]]
    opt_psisq1 <- mle_taupsi[["psisq"]]
    #
    not_rand_rows <- setdiff(1:nrow(Y), rand_rows)
    mle_taupsi <- BTF::find_signal_noise_mle(Y[not_rand_rows, , drop = F], U[not_rand_rows, ], V, R[not_rand_rows, , drop = F], sigmasq_tild)
    opt_tausq2 <- mle_taupsi[["tausq"]]
    opt_psisq2 <- mle_taupsi[["psisq"]]
    #
    rand_index <- c(sapply(rand_rows, function(rownum){(1:ncol(Y) - 1)*nrow(Y) + rownum}))
    not_rand_index <- c(sapply(not_rand_rows, function(rownum){(1:ncol(Y) - 1)*nrow(Y) + rownum}))
    perm_indx <- order(c(rand_index, not_rand_index)) - 1
  } else {
    #
    rand_cols <- sample(1:ncol(Y), round(ncol(Y)/2))
    mle_taupsi <- BTF::find_signal_noise_mle(Y[, rand_cols, drop = F], U, V[rand_cols, ], R[, rand_cols, drop = F], sigmasq_tild)
    opt_tausq1 <- mle_taupsi[["tausq"]]
    opt_psisq1 <- mle_taupsi[["psisq"]]
    #
    not_rand_cols <- setdiff(1:ncol(Y), rand_cols)
    mle_taupsi <- BTF::find_signal_noise_mle(Y[, not_rand_cols, drop = F], U, V[not_rand_cols, ], R[, not_rand_cols, drop = F], sigmasq_tild)
    opt_tausq2 <- mle_taupsi[["tausq"]]
    opt_psisq2 <- mle_taupsi[["psisq"]]
    #
    rand_index <- c(sapply(rand_cols, function(colnum){(colnum - 1)*nrow(Y) + 1:nrow(Y)}))
    not_rand_index <- c(sapply(not_rand_cols, function(colnum){(colnum - 1)*nrow(Y) + 1:nrow(Y)}))
    perm_indx <- order(c(rand_index, not_rand_index)) - 1
  }
  #
  opt_tausq <- c(opt_tausq1, opt_tausq2)
  opt_psisq <- c(opt_psisq1, opt_psisq2)
  #
  if(!is.null(Y1)){
    vR <- 1/((1/vecR) + (1/c(R1)))
  } else {
    vR <- vecR
  }
  linking_estimators <- BTF::rcpp_fabp_lin_reg(vecY, sigmasq_tild, opt_tausq, opt_psisq, vR, U, V, PermIndx = perm_indx)
  # Extract linking model estimators
  theta_hat <- linking_estimators[[2]]
  tau_hat <- linking_estimators[[1]]
  # Compute LOOCV and LOOR^2
  loocv <- mean((vecY - theta_hat)^2)
  loor2 <- cor(vecY, theta_hat)^2
  # Compute test t-statistics and corresponding UMP, FAB p-values
  if(is.null(Y1)){
    if(pool_sampling_var){
      tstat <- vecY/(sqrt(sigmasq_hat/vR))
    } else {
      tstat <- vecY/(sqrt(vecS^2/vR))
    }
  } else {
    if(pool_sampling_var){
      tstat <- vecY/sqrt(sigmasq_hat/vR)
    } else {
      vecSsq <- ((vecR-1)*(vecS^2) + c(R1-1)*(c(S1)^2))/(vecR + c(R1) - 2)
      tstat <- vecY/sqrt(vecSsq/vR)
    }
  }
  b <- 2*theta_hat*sqrt(sigmasq_tild/vR)/tau_hat
  if(is.null(Y1)){
    if(pool_sampling_var){
      p_values <- (1 - abs(pt(tstat, df = hat_nu) - pt(-tstat, df = hat_nu)))
      fabp_values <- (1 - abs(pt(tstat + b, df = hat_nu) - pt(-tstat, df = hat_nu)))
    } else {
      p_values <- (1 - abs(pt(tstat, df = vecR - 1) - pt(-tstat, df = vecR - 1)))
      fabp_values <- (1 - abs(pt(tstat + b, df = vecR - 1) - pt(-tstat, df = vecR - 1)))
    }
  } else {
    if(pool_sampling_var){
      p_values <- (1 - abs(pt(tstat, df = hat_nu) - pt(-tstat, df = hat_nu)))
      fabp_values <- (1 - abs(pt(tstat + b, df = hat_nu) - pt(-tstat, df = hat_nu)))
    } else {
      p_values <- (1 - abs(pt(tstat, df = vecR + c(R1) - 2) - pt(-tstat, df = vecR + c(R1) - 2)))
      fabp_values <- (1 - abs(pt(tstat + b, df = vecR + c(R1) - 2) - pt(-tstat, df = vecR + c(R1) - 2)))
    }
  }
  #
  return(data.frame(row = rep(row.names(Y), m), column = rep(colnames(Y), each = n),
                    observed = vecY, predicted = theta_hat, model_error = mean(linking_estimators[[3]]), theta_var = tau_hat, error_var = mean(c(sigmasq_hat, sigmasq_tild)),
                    statistic = tstat, guess = b,
                    p = p_values, fabp = fabp_values,
                    fdr_p = p.adjust(p_values, method = "BH"), fdr_fabp = p.adjust(fabp_values, method = "BH"),
                    mse = loocv, r2 = loor2))
}
#
