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
matrix_type <- c(0, 1, 0, 0, 2)
#
TensorList <- BTF::simulate_tensor_data(N, M, D, K, MatrixType = matrix_type, stn = 2)
#
sigmasq <- 1
R_true <- matrix(rnorm(D^2, sd = 1), ncol = D)
#
n <- 10
m <- 25
X <- kronecker(TensorList$V[1:m, 1:D], TensorList$U[1:n, 1:D])
eigs <- svd(X)[["u"]]
plist <- list()
for(tau in c(0.05, 0.25, 0.5, 0.75, 0.95)){
  L <- chol(tau*diag(rep(sigmasq, nrow(X))) + (1-tau)*eigs%*%t(eigs)) %>% t()
  #
  p_vals <- list()
  #
  for(t in 1:1000){
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
    pval_df <- BTF::fabp_lin_reg(Y = Y, S = S, R = matrix(r, n, m), U = TensorList$U[1:n, ], V = TensorList$V[1:m, ], pool_sampling_var = F)
    #
    p_vals[[paste0("", t)]] <- pval_df %>% dplyr::mutate(n = n*m, trial = t)
  }
  #
  p_vals <- data.table::rbindlist(p_vals)
  #
  wes_pal <- wesanderson::wes_palette("Zissou1", 5)
  #
  eval_quantiles <- seq(0, 1, length.out = 11)
  p1 <- data.frame(qp = c(quantile(p_vals$p, eval_quantiles), quantile(p_vals$fabp, eval_quantiles)),
                   type = c(rep("UMPU", 11), rep("FAB", 11)),
                   x = eval_quantiles) %>%
    dplyr::mutate(type = factor(type, levels = c("UMPU", "FAB"))) %>%
    ggplot2::ggplot() +
    geom_bar(aes(x = as.factor(round(x, 2)), y = qp, fill = type), stat = "identity", positio = "dodge") +
    scale_fill_manual(values = c(wes_pal[1], wes_pal[4]), name = "Test") +
    geom_point(aes(x = as.factor(round(x, 2)), y = x), colour = wes_pal[5], shape = 4) +
    theme_light() +
    labs(y = "Empirical quantile", x = "Theoretical quantile", title = paste0("t = ", tau)) +
    theme(legend.position = c(0.2, 0.7))
  #
  p3 <- p_vals %>%
    plyr::ddply(.(trial), function(x){
      thresh <- eval_quantiles[2:length(eval_quantiles)] / 5
      data.frame(UMPU = sapply(thresh, function(t){any(x$fdr_p < t)}),
                 FAB = sapply(thresh, function(t){any(x$fdr_fabp < t)}),
                 fdr = thresh)
    }) %>%
    reshape2::melt(id.vars = c("trial", "fdr"), variable.name = "type") %>%
    dplyr::group_by(type, fdr) %>%
    dplyr::summarise(mc_mean = mean(value), mc_err = sd(value) / sqrt(n())) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(type = factor(type, levels = c("UMPU", "FAB"))) %>%
    ggplot2::ggplot() +
    geom_bar(aes(x = as.factor(fdr), y = mc_mean, fill = type, group = type), stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(x = as.factor(fdr), ymin = mc_mean - mc_err, ymax = mc_mean + mc_err, colour = type, group = type), position = position_dodge(), colour = "black", size = 0.25) +
    scale_fill_manual(values = c(wes_pal[1], wes_pal[4]), name= "Test") +
    geom_point(aes(x = as.factor(fdr), y = fdr), colour = wes_pal[5], shape = 4) +
    theme_light() +
    labs(x = "Theoretical FDR", y = "Empirical FDR", title = "") +
    guides(fill = F)
  #
  plist[[paste0(tau, "")]] <- cowplot::plot_grid(p1, p3, nrow = 1)
}
cowplot::plot_grid(plotlist = plist, nrow = 5)
ggplot2::ggsave("../../figs/confounded_null_simulation_study_plot.pdf", width = 8, height = 15)
