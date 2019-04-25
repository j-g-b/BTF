#'
#'
#'
run_btf <- function(n_samples, tensor_list, embedding_dim = 20, inits = NULL, save_dir = NULL, thin = 10, burn = 10){
  #
  require(plyr)
  require(magrittr)
  #
  tensor_list %<>% magrittr::set_names(paste0("matrix", 0:(length(tensor_list) - 1)))
  matrix_type <- tensor_list %>%
                  plyr::laply(function(m){
                    if(all(m[!is.na(m)] %in% c(0, 1), na.rm = T)){
                      1
                    } else {
                      0
                    }
                  })
  #
  N <- nrow(tensor_list[["matrix0"]])
  M <- ncol(tensor_list[["matrix0"]])
  D <- embedding_dim
  K <- length(tensor_list)
  S <- n_samples
  #
  if(is.null(inits)){
    U0 <- matrix(0, nrow = N, ncol = D)
    V0 <- matrix(0, nrow = M, ncol = D)
    R0 <- plyr::alply(1:K, 1, function(k){
      rep(1, D) %>% diag()
    }) %>%
      magrittr::set_names(names(tensor_list))
    SS0 <- rep(0.5, K)
    LambdaU0 <- rep(1, D) %>% diag()
    LambdaV0 <- rep(1, D) %>% diag()
    MuU0 <- rep(0, D)
    MuV0 <- rep(0, D)
    Xi0 <- rep(0, D*D)
    Psi0 <- rep(1, D*D) %>% diag()
    IterStart <- 0
  } else {
    U0 <- inits[["U0"]]
    V0 <- inits[["V0"]]
    R0 <- inits[["R0"]]
    SS0 <- inits[["SS0"]]
    LambdaU0 <- inits[["LambdaU0"]]
    LambdaV0 <- inits[["LambdaV0"]]
    MuU0 <- inits[["MuU0"]]
    MuV0 <- inits[["MuV0"]]
    Xi0 <- inits[["Xi0"]]
    Psi0 <- inits[["Psi0"]]
    IterStart <- inits[["IterStart"]]
    burn <- 0
  }
  #
  if(is.null(save_dir)){
    save_dir <- "0"
  } else {
    if(!dir.exists(save_dir)){
      dir.create(save_dir)
    }
  }
  #
  btf_res <- BTF::BTF(TensorList = tensor_list, MatrixType = matrix_type,
                      U0 = U0, V0 = V0, R0 = R0, SS0 = SS0,
                      LambdaU0 = LambdaU0, LambdaV0 = LambdaV0,
                      MuU0 = MuU0, MuV0 = MuV0, Xi0 = Xi0, Psi0 = Psi0,
                      S = S, SaveDir = save_dir,
                      Thin = thin, Burn = burn, Start = IterStart)
  #
  btf_res[["U"]] %<>% magrittr::set_rownames(row.names(tensor_list[["matrix0"]]))
  btf_res[["V"]] %<>% magrittr::set_rownames(colnames(tensor_list[["matrix0"]]))
  #
  return(btf_res)

}