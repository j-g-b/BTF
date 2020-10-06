#'
#'
#'
run_btf <- function(n_samples, tensor_list, features_per_mod, inits = NULL, save_dir = NULL, thin = 10, burn = 10, matrix_type = NULL, method = "gibbs"){
  #
  require(plyr)
  require(magrittr)
  #
  tensor_list %<>% magrittr::set_names(paste0("matrix", 0:(length(tensor_list) - 1)))
  if(is.null(matrix_type)){
    matrix_type <- tensor_list %>%
      plyr::laply(function(m){
        if(all(m[!is.na(m)] %in% c(0, 1), na.rm = T)){
          1
        } else {
          0
        }
      })
  }
  #
  N <- nrow(tensor_list[["matrix0"]])
  M <- ncol(tensor_list[["matrix0"]])
  K <- length(tensor_list)
  d1 <- features_per_mod*K
  d2 <- features_per_mod*K
  S <- n_samples
  #
  if(is.null(inits)){
    interp_tensor_list <- lapply(1:length(tensor_list), function(i){
      res <- tensor_list[[i]] %>%
        apply(2, function(x){
          mx <- mean(x, na.rm = T)
          ifelse(is.na(x), mx, x)
        }) %>%
        apply(1, function(x){
          mx <- mean(x, na.rm = T)
          ifelse(is.na(x), mx, x)
        }) %>%
        t()
      if(matrix_type[i] == 1){
        res %>%
          pmax(0.025) %>%
          pmin(0.975) %>%
          qnorm()
      } else {
        res
      }
    })
    mat_list <- plyr::alply(1:K, 1, function(k){
      cat(paste0("\rDecomposing matrix ", k))
      U0 <- svd(interp_tensor_list[[k]], nu = features_per_mod)[["u"]]
      V0 <- svd(t(interp_tensor_list[[k]]), nu = features_per_mod)[["u"]]
      R0 <- t(U0)%*%(interp_tensor_list[[k]])%*%V0 + matrix(rnorm(features_per_mod^2, sd = 0.01), nrow = features_per_mod)
      list(U = U0, V = V0, R = R0)
    })
    U0 <- plyr::llply(mat_list, function(l){l[["U"]]}) %>% do.call(cbind, .)
    V0 <- plyr::llply(mat_list, function(l){l[["V"]]}) %>% do.call(cbind, .)
    R0 <- plyr::alply(1:K, 1, function(k){
            R <- mat_list[[k]][["R"]]
            zmat <- matrix(0, nrow = d1, ncol = d2)
            start <- (k-1)*features_per_mod + 1
            end <- k*features_per_mod
            zmat[start:end, start:end] <- R
            zmat
          }) %>%
          magrittr::set_names(names(tensor_list))
    SS0 <- sapply(1:length(tensor_list), function(i){0.0025})
    Mu0 <- sapply(1:length(tensor_list), function(i){0})
    IterStart <- 0
  } else {
    U0 <- inits[["U0"]]
    V0 <- inits[["V0"]]
    R0 <- inits[["R0"]]
    SS0 <- inits[["SS0"]]
    Mu0 <- inits[["Mu0"]]
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
  if(method == "gibbs"){
    btf_res <- BTF::BTF(TensorList = tensor_list, MatrixType = matrix_type,
                        U0 = U0, V0 = V0, R0 = R0, SS0 = SS0, Mu0 = Mu0,
                        S = S, SaveDir = save_dir,
                        Thin = thin, Burn = burn, Start = IterStart)
  } else if(method == "em"){
    btf_res <- BTF::EM(TensorList = tensor_list, MatrixType = matrix_type,
                        U0 = U0, V0 = V0, R0 = R0, SS0 = SS0, Mu0 = Mu0,
                        S = S)
  }
  #
  btf_res[["U"]] %<>% magrittr::set_rownames(row.names(tensor_list[["matrix0"]]))
  btf_res[["V"]] %<>% magrittr::set_rownames(colnames(tensor_list[["matrix0"]]))
  #
  return(btf_res)

}
