#
library(tidyverse)
library(magrittr)
library(plyr)
#
N <- 100
M <- 5000
r <- 5
D <- 5
K <- 5
n <- 10
m <- 25
matrix_type <- c(0, 1, 0, 0, 2)
#
TensorList <- BTF::simulate_tensor_data(N, M, D, K, MatrixType = matrix_type, stn = 2)
#
sigmasq <- 1
R_true <- matrix(rnorm(D^2, sd = 1), ncol = D)
design_eigs <- svd(kronecker(TensorList$V[1:m, 1:D], TensorList$U[1:n, 1:D]))[["u"]][, 1:D^2]
#
fdr_list <- list()
for(j in 1:500){
  #
  X <- matrix(rnorm(m*D*D*n), nrow = m*n)
  eigs <- svd(X)[["u"]][, 1:D^2]
  tau <- 0.05
  L <- chol(tau*diag(rep(sigmasq, nrow(X))) + (1-tau)*eigs%*%t(eigs)) %>% t()
  tv <- sum((t(eigs)%*%design_eigs)^2)
  #
  p_vals <- list()
  #
  for(t in 1:100){
    #
    print(t)
    #
    theta_true <- rep(0, n*m)
    #
    Y <- array(dim = c(n, m, r))
    for(i in 1:r){
      Y[, , i] <- L%*%rnorm(length(theta_true), mean = theta_true, sd = 1)
    }
    S <- apply(Y, c(1, 2), sd)
    Y <- apply(Y, c(1, 2), mean) %>%
      magrittr::set_rownames(paste0("", 1:n)) %>%
      magrittr::set_colnames(paste0("", 1:m))
    #
    pval_df <- BTF::fabp_lin_reg(Y = Y, S = S, R = matrix(r, n, m), U = TensorList$U[1:n, ], V = TensorList$V[1:m, ])
    #
    p_vals[[paste0("", t)]] <- pval_df %>% dplyr::mutate(n = n*m, trial = t)
  }
  #
  fdr_list[[j]] <- data.table::rbindlist(p_vals) %>%
    plyr::ddply(.(trial), function(x){
      thresh <- 0.1
      data.frame(UMPU = sapply(thresh, function(th){any(x$fdr_p < th)}),
                 FAB = sapply(thresh, function(th){any(x$fdr_fabp < th)}),
                 fdr = thresh)
    }) %>%
    reshape2::melt(id.vars = c("trial", "fdr"), variable.name = "type") %>%
    dplyr::group_by(type, fdr) %>%
    dplyr::summarise(mc_mean = mean(value), mc_err = sd(value) / sqrt(n())) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(trial = j, tv = tv)
}
#
fdr_df <- data.table::rbindlist(fdr_list)
#
wes_pal <- wesanderson::wes_palette("Zissou1", 5)
#
fdr_df %>%
  reshape2::dcast(trial + tv~type, value.var = "mc_mean") %>%
  dplyr::mutate(fab_mean = mean(FAB), umpu_mean = mean(UMPU)) %>%
  ggplot2::ggplot() +
  geom_hex(aes(x = UMPU, y = FAB), bins = 50) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(x = 0.1, y = 0.1, shape = 4, size = 10, colour = wes_pal[5]) +
  geom_hline(aes(yintercept = fab_mean), linetype = 1, colour = wes_pal[4]) +
  geom_vline(aes(xintercept = umpu_mean), linetype = 1, colour = wes_pal[1]) +
  scale_fill_distiller(palette = "Greys") +
  theme_light() +
  guides(fill = F) +
  labs(x = "UMPU empirical FDR", y = "FAB empirical FDR") +
  coord_equal()
#
ggplot2::ggsave("../../figs/model_misspec_null.pdf", height = 6, width = 6)

