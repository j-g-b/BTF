#
library(tidyverse)
library(magrittr)
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
#
p_vals <- list()
#
for(t in 1:10000){
  #
  print(t)
  #
  theta_true <- rep(0, n*m)
  #
  Y <- array(dim = c(n, m, r))
  for(i in 1:r){
    Y[, , i] <- rnorm(length(theta_true), mean = theta_true, sd = sqrt(sigmasq))
  }
  S <- apply(Y, c(1, 2), sd)
  Y <- apply(Y, c(1, 2), mean) %>% 
    magrittr::set_rownames(paste0("", 1:n)) %>% 
    magrittr::set_colnames(paste0("", 1:m))
  #
  pval_df <- BTF::fabp_lin_reg(Y = Y, S = S, R = matrix(r, n, m), U = svd(TensorList$U)[["u"]][1:n, 1:min(n, D)], V = svd(TensorList$V)[["u"]][1:m, 1:min(m, D)], snr = 2)
  #
  p_vals[[paste0("", t)]] <- pval_df %>% dplyr::mutate(n = n*m, trial = t)
}
#
p_vals <- data.table::rbindlist(p_vals)
#
wes_pal <- wesanderson::wes_palette("Zissou1", 5)
#
p1 <- data.frame(x = p_vals$p) %>%
  ggplot2::ggplot() +
  geom_histogram(aes(x = x, y = ..density..), breaks = seq(0, 1, length.out = 20), colour = "white", fill = wes_pal[2]) +
  theme_light() +
  labs(y = "Frequency UMPU", x = "")
p2 <- data.frame(x = p_vals$fabp) %>%
  ggplot2::ggplot() +
  geom_histogram(aes(x = x, y = ..density..), breaks = seq(0, 1, length.out = 20), colour = "white", fill = wes_pal[3]) +
  theme_light() +
  labs(y = "Frequency FAB", x = "")
#
p3 <- p_vals %>%
        dplyr::group_by(trial) %>%
        dplyr::summarise(UMPU = any(fdr_p < 0.1),
                         FAB = any(fdr_fabp < 0.1)) %>%
        dplyr::ungroup() %>%
        reshape2::melt(id.vars = "trial", variable.name = "type") %>%
        dplyr::group_by(type) %>%
        dplyr::summarise(mc_mean = mean(value), mc_err = sd(value) / sqrt(n())) %>%
        dplyr::ungroup() %>%
        ggplot2::ggplot() +
        geom_bar(aes(x = type, y = mc_mean, fill = type), stat = "identity") +
        geom_errorbar(aes(x = type, ymin = mc_mean - mc_err, ymax = mc_mean + mc_err), width = 0.2) +
        scale_fill_manual(values = c(wes_pal[2], wes_pal[3])) +
        theme_light() +
        geom_hline(yintercept = 0.1, size = 0.5, colour = wes_pal[5], alpha = 0.5) +
        labs(x = "", y = "FDR") +
        guides(fill = F)
#
cowplot::plot_grid(p1, p2, p3, nrow = 1)
ggplot2::ggsave("figs/null_simulation_study_plot.pdf", width = 10, height = 3.5)
