#
library(tidyverse)
library(magrittr)
#
N <- 100
M <- 5000
n <- 5
m <- 100
r <- 6
D <- 10
K <- 5
matrix_type <- c(0, 1, 0, 0, 2)
#
TensorList <- BTF::simulate_tensor_data(N, M, D, K, MatrixType = matrix_type, stn = 2)
#
sigmasq <- 1/rgamma(n*m, 10, 10)
psisq <- 1
R_true <- matrix(rnorm(D^2, sd = 1), ncol = D)
#
psiseq <- seq(0, 1, length.out = 6)
#
plist <- list()
plist2 <- list()
plist3 <- list()
#
var_component <- sqrt(psisq)*rnorm(n*m)
var_component <- var_component / sd(var_component)
linear_component <- TensorList$U[1:n, ]%*%R_true%*%t(TensorList$V[1:m, ])
linear_component <- linear_component / sd(linear_component)
linear_component <- c(linear_component)
#
wes_pal <- wesanderson::wes_palette("Zissou1", 5)
#
for(t in psiseq){
  #
  rank_p <- data.frame()
  rank_fabp <- data.frame()
  rank_osp <- data.frame()
  #
  for(j in 1:200){
    #
    print(j)
    #
    theta_true <- sqrt(t)*linear_component + sqrt((1 - t))*var_component
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
    pval_df <- BTF::fabp_lin_reg(Y = Y, S = S, R = matrix(r, n, m), U = TensorList$U[1:n, ], V = TensorList$V[1:m, ], pool_sampling_var = F)
    #
    rank_p <- rbind(rank_p, data.frame(p = pval_df$p, fdr = pval_df$fdr_p, rank = rank(pval_df$p), type = "UMPU", trial = j))
    rank_fabp <- rbind(rank_fabp, data.frame(p = pval_df$fabp, fdr = pval_df$fdr_fabp, rank = rank(pval_df$fabp), type = "FAB", trial = j))
    rank_osp <- rbind(rank_osp, data.frame(p = ifelse(theta_true < 0, pt(pval_df$statistic, df = (r - 1)), 1 - pt(pval_df$statistic, df = (r - 1)))) %>%
                                  dplyr::mutate(fdr = p.adjust(p, method = "BH"), rank = rank(p), type = "OS", trial = j))
    #
  }
  #
  if(abs(t) < 1e-6){
    plist2[[paste0(t, "")]] <- rbind(rank_p, rank_fabp, rank_osp)  %>%
      dplyr::mutate(type = factor(type, levels = c("UMPU", "FAB", "OS"))) %>%
      dplyr::filter(fdr < 0.1) %>%
      dplyr::mutate(fdr_bin = cut(fdr, breaks = c(-Inf, 0.001, 0.005, 0.01, 0.05, 0.1))) %>%
      dplyr::mutate(fdr_bin = gsub("]", "", gsub(".*,", "",fdr_bin))) %>%
      dplyr::group_by(type, fdr_bin) %>%
      dplyr::summarise(n_discov = n()) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(type) %>%
      dplyr::mutate(n_discov = cumsum(n_discov)) %>%
      dplyr::ungroup() %>%
      ggplot2::ggplot() +
      geom_bar(aes(x = fdr_bin, y = n_discov, fill = type), stat = "identity", position = "dodge") +
      theme_bw() +
      scale_fill_manual(name = "Test", labels = c("UMPU", "FAB", "OS"), values = c(wes_pal[1], wes_pal[4], "grey20")) +
      theme_light() +
      labs(title = bquote("" ~ tau^2 ~ " = " ~ .(round(psisq*(1 - t), 2))), x = "FDR", y = "# of discoveries") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1),
            legend.position = c(0.2, 0.7))
  } else {
    plist2[[paste0(t, "")]] <- rbind(rank_p, rank_fabp, rank_osp)  %>%
      dplyr::mutate(type = factor(type, levels = c("UMPU", "FAB", "OS"))) %>%
      dplyr::filter(fdr < 0.1) %>%
      dplyr::mutate(fdr_bin = cut(fdr, breaks = c(-Inf, 0.001, 0.005, 0.01, 0.05, 0.1))) %>%
      dplyr::mutate(fdr_bin = gsub("]", "", gsub(".*,", "",fdr_bin))) %>%
      dplyr::group_by(type, fdr_bin) %>%
      dplyr::summarise(n_discov = n()) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(type) %>%
      dplyr::mutate(n_discov = cumsum(n_discov)) %>%
      dplyr::ungroup() %>%
      ggplot2::ggplot() +
      geom_bar(aes(x = fdr_bin, y = n_discov, fill = type), stat = "identity", position = "dodge") +
      theme_bw() +
      scale_fill_manual(name = "Test", labels = c("UMPU", "FAB", "OS"), values = c(wes_pal[1], wes_pal[4], "grey20")) +
      theme_light() +
      labs(title = bquote("" ~ tau^2 ~ " = " ~ .(round(psisq*(1 - t), 2))), x = "FDR", y = "# of discoveries") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none")
  }
}
#
cowplot::plot_grid(plotlist = plist2, nrow = 2)
ggplot2::ggsave("../../figs/model_precision_sim_ndiscov_plot.pdf", height = 7, width = 10)
