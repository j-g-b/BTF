#'
#'
#'
align_samples <- function(sample_list){
  #
  require(plyr)
  require(magrittr)
  #
  S <- dim(sample_list[["U"]])[1]
  #
  U0 <- sample_list[["U"]][S, ,] %>% apply(2, function(u){ u / sum(u^2)})
  V0 <- sample_list[["V"]][S, ,] %>% apply(2, function(v){ v / sum(v^2)})
  #
  aligned_samples <- plyr::alply(1:S, 1, function(s){
    #
    U <- sample_list[["U"]][s,,]
    V <- sample_list[["V"]][s,,]
    R <- sample_list[["R"]][s,,,]
    #
    col_norms_U <- apply(U, 2, function(u){sqrt(sum(u^2))})
    col_norms_V <- apply(V, 2, function(v){sqrt(sum(v^2))})
    #
    U <- apply(U, 2, function(u){ u / sqrt(sum(u^2))})
    V <- apply(V, 2, function(v){ v / sqrt(sum(v^2))})
    svd_U <- svd(t(U)%*%U0)
    svd_V <- svd(t(V)%*%V0)
    #
    U_Q <- svd_U[["u"]]%*%t(svd_U[["v"]])
    V_Q <- svd_V[["u"]]%*%t(svd_V[["v"]])
    #
    U_align <- U%*%U_Q
    V_align <- V%*%V_Q
    R_align <- plyr::aaply(R, 1, function(r){
      t(U_Q)%*%diag(col_norms_U)%*%r%*%diag(col_norms_V)%*%t(V_Q)
    })
    #
    list(U = U_align, V = V_align, R = R_align)
  })
  #
  sample_list[["U"]] <- plyr::laply(aligned_samples, function(l){ l[["U"]] })
  sample_list[["V"]] <- plyr::laply(aligned_samples, function(l){ l[["V"]] })
  sample_list[["R"]] <- plyr::laply(aligned_samples, function(l){ l[["R"]] })
  #
  return(sample_list)
}
