#'
#'
#'
align_samples <- function(sample_list){
  #
  require(plyr)
  require(magrittr)
  #
  iter <- 5
  S <- dim(sample_list[["U"]])[1]
  #
  U0 <- sample_list[["U"]][S, ,] %>% apply(2, function(u){ u / sqrt(sum(u^2))})
  V0 <- sample_list[["V"]][S, ,] %>% apply(2, function(v){ v / sqrt(sum(v^2))})
  #
  for(i in 1:iter){
    #
    cat(paste0("\rAligning samples iter: ", i))
    #
    aligned_samples <- plyr::alply(1:S, 1, function(s){
      #
      U <- sample_list[["U"]][s,,]
      V <- sample_list[["V"]][s,,]
      R <- sample_list[["R"]][s,,,]
      #
      svd_U <- svd(t(U)%*%U0)
      svd_V <- svd(t(V)%*%V0)
      #
      U_Q <- svd_U[["u"]]%*%t(svd_U[["v"]])
      V_Q <- svd_V[["u"]]%*%t(svd_V[["v"]])
      #
      U_align <- U%*%U_Q
      V_align <- V%*%V_Q
      #
      col_norms_U <- apply(U_align, 2, function(u){sqrt(sum(u^2))})
      col_norms_V <- apply(V_align, 2, function(v){sqrt(sum(v^2))})
      #
      U_align <- apply(U_align, 2, function(u){ u / sqrt(sum(u^2))})
      V_align <- apply(V_align, 2, function(v){ v / sqrt(sum(v^2))})
      #
      R_align <- plyr::aaply(R, 1, function(r){
        diag(col_norms_U)%*%t(U_Q)%*%r%*%V_Q%*%diag(col_norms_V)
      })
      #
      list(U = U_align, V = V_align, R = R_align)
    })
    #
    sample_list[["U"]] <- plyr::laply(aligned_samples, function(l){ l[["U"]] })
    sample_list[["V"]] <- plyr::laply(aligned_samples, function(l){ l[["V"]] })
    sample_list[["R"]] <- plyr::laply(aligned_samples, function(l){ l[["R"]] })
    #
    U_prev <- U0
    #
    U0 <- apply(sample_list[["U"]], c(2, 3), mean) %>% apply(2, function(u){ u / sqrt(sum(u^2))})
    V0 <- apply(sample_list[["V"]], c(2, 3), mean) %>% apply(2, function(v){ v / sqrt(sum(v^2))})
    #
  }
  #
  return(sample_list)
}
