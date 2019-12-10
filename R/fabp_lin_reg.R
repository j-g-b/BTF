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
fabp_lin_reg <- function(Y, S, R, U, V){
  # Get experimental dimensions
  vecY <- c(Y)
  vecR <- c(R)
  vecS <- c(S)
  m <- nrow(V)
  n <- nrow(U)
  # Split data and estimate sigmasq_hat and sigmasq_tild
  hat_indices <- sample(m*n, round(m*n / 2))
  tild_indices <- setdiff(1:(m*n), hat_indices)
  hat_nu <- sum(vecR[hat_indices])
  tild_nu <- sum(vecR[tild_indices])
  sigmasq_hat <- sum((R[hat_indices] - 1)*S[hat_indices]^2) / hat_nu
  sigmasq_tild <- sum((R[tild_indices] - 1)*S[tild_indices]^2) / tild_nu
  #
  linking_estimators <- BTF::rcpp_fabp_lin_reg(vecY, sigmasq_hat, vecR, U, V)
  # Extract linking model estimators
  theta_hat <- linking_estimators[[2]]
  tau_hat <- linking_estimators[[1]]
  # Compute test t-statistics and corresponding UMP, FAB p-values
  tstat <- vecY/(sqrt(sigmasq_hat/vecR))
  b <- 2*theta_hat*sqrt(sigmasq_tild/vecR)/tau_hat
  #b[abs(b) > 4] <- sign(b)*4
  p_values <- (1 - abs(pt(tstat, df = hat_nu) - pt(-tstat, df = hat_nu))) %>% p.adjust(method = "BH")
  fabp_values <- (1 - abs(pt(tstat + b, df = hat_nu) - pt(-tstat, df = hat_nu))) %>% p.adjust(method = "BH")
  #
  return(data.frame(row = rep(row.names(Y), m), column = rep(colnames(Y), each = n), 
                    observed = vecY, predicted = theta_hat, statistic = tstat, guess = b,
                    p = p_values, fabp = fabp_values))
}