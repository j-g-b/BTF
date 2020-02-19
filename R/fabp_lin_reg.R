#
#' Compute FAB p-values under linear regression linking model
#'
#' @description Comutes FAB p-values based on a linear regression linking model for matrix data given row and column covariates.
#'
#' @param Y matrix of averages of experimental readout values over R replicates
#' @param S matrix of standard errors of experimental readout values over R replicates
#' @param R number of replicates per readout value (can be matrix or scalar; if scalar assumes that number of replicates was the same for all data values in Y)
#' @param X design matrix (the Kronecker product of cell line and gene covariates)
#'
#' @return A data.frame of FAB p-values and the standard UMP p-values, one for each entry in Y
#'
#' @export fabp_lin_reg
#'
fabp_lin_reg <- function(Y, S, R, U, V, contrasts = NULL, snr = 0.5){
  # Get experimental dimensions
  vecY <- c(Y)
  vecR <- c(R)
  m <- nrow(V)
  n <- nrow(U)
  p <- ncol(V)
  q <- ncol(U)
  # Split data and estimate sigmasq_hat and sigmasq_tild
  if(is.matrix(S)){
    hat_indices <- sample(m*n, round(m*n / 2))
    tild_indices <- setdiff(1:(m*n), hat_indices)
    hat_nu <- sum(vecR[hat_indices] - 1)
    tild_nu <- sum(vecR[tild_indices] - 1)
    sigmasq_hat <- sum((R[hat_indices] - 1)*S[hat_indices]^2) / hat_nu
    sigmasq_tild <- sum((R[tild_indices] - 1)*S[tild_indices]^2) / tild_nu
  } else if(is.null(dim(S))){
    sigmasq_hat <- S^2
    sigmasq_tild <- S^2
    hat_nu <- sum(vecR)
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
  linking_estimators <- BTF::rcpp_fabp_lin_reg(vecY, sigmasq_hat, opt_tausq, opt_psisq, vecR, U, V, PermIndx = perm_indx)
  # Extract linking model estimators
  theta_hat <- linking_estimators[[2]]
  tau_hat <- linking_estimators[[1]]
  # Compute LOOCV and LOOR^2
  loocv <- mean((vecY - theta_hat)^2)
  loor2 <- cor(vecY, theta_hat)^2
  # Compute test t-statistics and corresponding UMP, FAB p-values
  tstat <- vecY/(sqrt(sigmasq_hat/vecR))
  b <- 2*theta_hat*sqrt(sigmasq_tild/vecR)/tau_hat
  #b[abs(b) > 4] <- sign(b)*4
  if(!is.null(dim(S))){
    p_values <- (1 - abs(pt(tstat, df = hat_nu) - pt(-tstat, df = hat_nu)))
    fabp_values <- (1 - abs(pt(tstat + b, df = hat_nu) - pt(-tstat, df = hat_nu)))
  } else {
    p_values <- (1 - abs(pnorm(tstat) - pnorm(-tstat)))
    fabp_values <- (1 - abs(pnorm(tstat + b) - pnorm(-tstat)))
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
