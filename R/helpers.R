#
simulate_tensor_data <- function(N, M, D, K, MatrixType, rho = 0){
  #
  SigmaSq <- 1/rgamma(K, 5/2, 1/2)
  #
  U <- rnorm(N*D) %>%
    matrix(nrow = N) %>%
    apply(1, function(u){u / sqrt(sum(u^2))}) %>%
    t()
  V <- rnorm(M*D) %>%
    matrix(nrow = M) %>%
    apply(1, function(v){v / sqrt(sum(v^2))}) %>%
    t()
  R <- plyr::aaply(1:(D*D), 1, function(d){
    W <- matrix(rho, nrow = K, ncol = K) + diag(1 - rho, nrow = K)
    R <- (t(chol(W))%*%rnorm(K))
  }) %>%
    plyr::alply(2, function(x){matrix(x, nrow = D, ncol = D)})
  #
  X <- plyr::llply(1:K, function(k){
    if(MatrixType[k] == 0){
      (U%*%R[[k]]%*%t(V)) %>%
        magrittr::add(matrix(rnorm(N*M, sd = sqrt(SigmaSq[k])), nrow = N)) %>%
        magrittr::set_rownames(paste0("CL", 1:N)) %>%
        magrittr::set_colnames(paste0("GN", 1:M))
    } else if(MatrixType[k] == 1){
      (U%*%R[[k]]%*%t(V)) %>%
        magrittr::add(matrix(rnorm(N*M, sd = 1), nrow = N)) %>%
        apply(2, function(Z){ifelse(rbernoulli(length(Z), pnorm(Z)), 1, 0)}) %>%
        magrittr::set_rownames(paste0("CL", 1:N)) %>%
        magrittr::set_colnames(paste0("GN", 1:M))
    } else if(MatrixType[k] == 2){
      (U%*%R[[k]]%*%t(V)) %>%
        magrittr::add(matrix(rnorm(N*M, sd = sqrt(SigmaSq[k])), nrow = N)) %>%
        apply(2, function(Z){ifelse(Z < 0, 0, Z)}) %>%
        magrittr::set_rownames(paste0("CL", 1:N)) %>%
        magrittr::set_colnames(paste0("GN", 1:M))
    }
  })
  #
  return(list(U = U, V = V, R = R, SigmaSq = SigmaSq, X = X))
}
