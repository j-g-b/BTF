#
samps_dir <- "../../alt_samps"
#
library(tidyverse)
library(magrittr)
#
samps <- BTF::collect_samples(samps_dir, last = 100) %>% BTF::align_samples()
#
ras_mut <- readr::read_csv("ras_mut_status.csv") %>%
  reshape2::melt(id.vars = "status", variable.name = "cl") %>%
  dplyr::group_by(cl) %>%
  dplyr::filter(value[status == "is AML"] == 1) %>%
  dplyr::ungroup() %>%
  dplyr::filter(status != "is AML") %>%
  dplyr::group_by(cl) %>%
  dplyr::summarise(is_ras_mut = any(value == 1)) %>%
  dplyr::ungroup()
#
n_perm <- 2000
p_val_list <- list()
for(i in 1:n_perm){
  ras_mut %<>% 
    dplyr::group_by(is_ras_mut) %>% 
    dplyr::mutate(is_ras_mut_perm = ifelse(sample(1:n(), n())%%2 == 0, T, F)) %>%
    dplyr::ungroup()
  table <- readr::read_csv("focused_screen.csv") %>%
    reshape2::melt(id.vars = "Gene", variable.name = "cl") %>%
    dplyr::right_join(ras_mut) %>%
    dplyr::group_by(Gene, is_ras_mut_perm) %>%
    dplyr::summarise(r = sum(!is.na(value)),
                     s = sd(value, na.rm = T),
                     m = mean(value, na.rm = T)) %>%
    dplyr::ungroup() %>%
    plyr::ddply(.(Gene), function(d){
      data.frame(r = d$r[d$is_ras_mut_perm],
                 s = d$s[d$is_ras_mut_perm],
                 m = d$m[d$is_ras_mut_perm],
                 r1 = d$r[!d$is_ras_mut_perm],
                 s1 = d$s[!d$is_ras_mut_perm],
                 m1 = d$m[!d$is_ras_mut_perm])
    })
  #
  Ybar <- table$m %>%
    matrix(nrow = 1) %>%
    magrittr::set_colnames(table$Gene)
  #
  Ybar1 <- table$m1 %>%
    matrix(nrow = 1) %>%
    magrittr::set_colnames(table$Gene)
  #
  S <- table$s %>%
    matrix(nrow = 1) %>%
    magrittr::set_colnames(table$Gene)
  #
  S1 <- table$s1 %>%
    matrix(nrow = 1) %>%
    magrittr::set_colnames(table$Gene)
  #
  R <- table$r %>%
    matrix(nrow = 1) %>%
    magrittr::set_colnames(table$Gene)
  #
  R1 <- table$r1 %>%
    matrix(nrow = 1) %>%
    magrittr::set_colnames(table$Gene)
  #
  if(i == 1){
    gene_map <- readr::read_csv("gene_map.csv")
    V <- apply(samps$V, c(2,3), mean) %>% magrittr::set_rownames(gene_map[["gene"]])
    #
    gns_to_use <- intersect(row.names(V), colnames(S))
    V <- V[gns_to_use, ]
  }
  #
  Ybar <- Ybar[, gns_to_use, drop = F] %>% magrittr::set_rownames("diff_score")
  S <- S[, gns_to_use, drop = F] %>% magrittr::set_rownames("diff_score")
  R <- R[, gns_to_use, drop = F] %>% magrittr::set_rownames("diff_score")
  #
  Ybar1 <- Ybar1[, gns_to_use, drop = F] %>% magrittr::set_rownames("diff_score")
  S1 <- S1[, gns_to_use, drop = F] %>% magrittr::set_rownames("diff_score")
  R1 <- R1[, gns_to_use, drop = F] %>% magrittr::set_rownames("diff_score")
  #
  p_vals <- BTF::fabp_lin_reg(Y = Ybar, S = S, R = R, 
                              U = matrix(1, nrow = 1, ncol = 1), V = V, 
                              pool_sampling_var = F,
                              Y1 = Ybar1, S1 = S1, R1 = R1)
  #
  p_val_list[[i]] <- p_vals %>% dplyr::mutate(perm = i)
}
#
p_val_df <- p_val_list %>% do.call(rbind, .)
#
wes_pal <- wesanderson::wes_palette("Zissou1", 5)
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
ggplot2::ggsave("../../figs/wang_focused_null_example.pdf", height = 5, width = 10)

