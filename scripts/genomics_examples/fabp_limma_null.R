#
samps_dir <- "../../alt_samps"
#
library(limma)
library(affy)
library(tidyverse)
library(magrittr)
#
targets <- readTargets()
x <- read.ilmn(files="probe profile.txt", ctrlfiles="control probe profile.txt", other.columns="Detection")
y <- neqc(x)
expressed <- rowSums(y$other$Detection < 0.05) >= 3
y <- y[expressed,]
ct <- factor(targets$CellType)
#
samps <- BTF::collect_samples(samps_dir, last = 100) %>% BTF::align_samples()
#
wes_pal <- wesanderson::wes_palette("Zissou1", 5)
#
contrasts <- matrix(c("ML", "LP", "ML", "MS", "MS", "LP"), nrow = 3)
y$E <- y$E[, ct != "Stroma"]
ct <- ct[ct != "Stroma"]
#
n_perm <- 500
#
plot_list <- list()
p_val_list <- list()
#
for(i in 1:n_perm){
  cat(paste0("\r", i))
  trt1 <- y$E[, ct == "ML"]
  trt2 <- y$E[, ct == "LP"]
  trt3 <- y$E[, ct == "MS"]
  rsamp <- cbind(sample(1:3, 2), sample(1:3, 2), sample(1:3, 2))
  trt_perm <- cbind(trt1[, rsamp[1, 1]], trt2[, rsamp[1, 2]], trt3[, rsamp[1, 3]])
  ctrl_perm <- cbind(trt1[, rsamp[2, 1]], trt2[, rsamp[2, 2]], trt3[, rsamp[2, 3]])
  #
  Y <- rowMeans(trt_perm) %>%
    matrix(ncol = 1) %>%
    magrittr::set_rownames(y$genes$SYMBOL) %>%
    magrittr::set_colnames("diff_score")
  #
  Y1 <- rowMeans(ctrl_perm) %>%
    matrix(ncol = 1) %>%
    magrittr::set_rownames(y$genes$SYMBOL) %>%
    magrittr::set_colnames("diff_score")
  #
  S <- apply(trt_perm, 1, sd) %>%
    matrix(ncol = 1) %>%
    magrittr::set_rownames(y$genes$SYMBOL) %>%
    magrittr::set_colnames("diff_score")
  #
  S1 <- apply(ctrl_perm, 1, sd) %>%
    matrix(ncol = 1) %>%
    magrittr::set_rownames(y$genes$SYMBOL) %>%
    magrittr::set_colnames("diff_score")
  #
  R <- rep(3, nrow(Y)) %>% 
    matrix(ncol = 1) %>%
    magrittr::set_rownames(y$genes$SYMBOL) %>%
    magrittr::set_colnames("diff_score")
  #
  R1 <- rep(3, nrow(Y)) %>% 
    matrix(ncol = 1) %>%
    magrittr::set_rownames(y$genes$SYMBOL) %>%
    magrittr::set_colnames("diff_score")
  #
  if(i == 1){
    gene_map <- readr::read_csv("gene_map.csv")
    V <- apply(samps$V, c(2,3), mean) %>% magrittr::set_rownames(gene_map[["gene"]])
    gns_to_use <- intersect(row.names(V), row.names(S))
    V <- V[gns_to_use, ]
  }
  #
  p_vals <- BTF::fabp_lin_reg(Y = t(Y[gns_to_use, , drop=F]), 
                              S = t(S[gns_to_use, , drop=F]), 
                              R = t(R[gns_to_use, , drop=F]), 
                              V = V, U = matrix(1, nrow = 1, ncol = 1), 
                              pool_sampling_var = T,
                              Y1 = t(Y1[gns_to_use, , drop=F]), 
                              S1 = t(S1[gns_to_use, , drop=F]), 
                              R1 = t(R1[gns_to_use, , drop=F]))
  #
  rank_p <- data.frame(p = p_vals$p, fdr = p_vals$fdr_p, rank = rank(p_vals$p), type = "UMPU", observed = p_vals$statistic)
  rank_fabp <- data.frame(p = p_vals$fabp, fdr = p_vals$fdr_fabp, rank = rank(p_vals$fabp), type = "FAB", observed = p_vals$statistic)
  #
  p_val_list[[i]] <- p_vals %>% dplyr::mutate(perm = i)
}
#
p_val_df <- p_val_list %>% do.call(rbind, .)
#
p1 <- p_val_df %>%
  plyr::ddply(.(perm), function(df){
    data.frame(x = rep(seq(0, 1, length.out = 100), 3),
               p = c(ecdf(df$p)(seq(0, 1, length.out = 100)),
                     ecdf(df$fabp)(seq(0, 1, length.out = 100)),
                     ecdf(runif(nrow(df)))(seq(0, 1, length.out = 100))),
               type = rep(c("UMPU", "FAB", "Uniform"), each = 100))
  }) %>%
  dplyr::group_by(x, type) %>%
  dplyr::summarise(m = mean(p), lwr = quantile(p, 0.025), upr = quantile(p, 0.975)) %>%
  dplyr::ungroup() %>%
  ggplot2::ggplot() +
  geom_ribbon(aes(x = x, ymin = lwr, ymax = upr, fill = type), alpha = 0.2) +
  geom_line(aes(x = x, y = m, colour = type), size = 1) +
  theme_light() +
  scale_colour_manual(values = c(wes_pal[4], wes_pal[2], wes_pal[5]), name = "") +
  scale_fill_manual(values = c(wes_pal[4], wes_pal[2], wes_pal[5]), name = "") +
  labs(x = "p", y = "Empirical percentile") +
  theme(legend.position = c(0.2, 0.7))
#
p2 <- p_val_df %>%
  dplyr::group_by(perm) %>%
  dplyr::summarise(UMPU = ifelse(any(fdr_p < 0.1), 1, 0),
                   FAB = ifelse(any(fdr_fabp < 0.1), 1, 0)) %>%
  dplyr::ungroup() %>%
  reshape2::melt(id.vars = "perm", variable.name = "type") %>%
  dplyr::group_by(type) %>%
  dplyr::summarise(mc_mean = mean(value),
                   mc_sd = sd(value) / sqrt(n())) %>%
  ggplot2::ggplot() +
  geom_bar(aes(x = type, y = mc_mean, fill = type), stat = "identity") +
  geom_errorbar(aes(x = type, ymin = mc_mean - mc_sd, ymax = mc_mean + mc_sd), width = 0.2) +
  geom_hline(yintercept = 0.1, colour = wes_pal[5]) +
  scale_fill_manual(values = c(wes_pal[2], wes_pal[4]), name = "Type") +
  theme_light() +
  guides(fill = F) +
  labs(x = "Test", y = "Empirical FDR")
#
cowplot::plot_grid(p1, p2)
#
ggplot2::ggsave("../../figs/limma_null_example.pdf", width = 10, height = 5)      

