#
samps_dir <- "../../alt_samps"
#
library(tidyverse)
library(magrittr)
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
table <- readr::read_csv("focused_screen.csv") %>%
  reshape2::melt(id.vars = "Gene", variable.name = "cl") %>%
  dplyr::right_join(ras_mut) %>%
  dplyr::group_by(Gene, is_ras_mut) %>%
  dplyr::summarise(r = sum(!is.na(value)),
                   s = sd(value, na.rm = T),
                   m = mean(value, na.rm = T)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(r = mean(r),
                   s = sqrt(sum(s^2)),
                   m = m[is_ras_mut] - m[!is_ras_mut]) %>%
  dplyr::ungroup()
#
Ybar <- table$m %>%
  matrix(nrow = 1) %>%
  magrittr::set_colnames(table$Gene)
#
S <- table$s %>%
  matrix(nrow = 1) %>%
  magrittr::set_colnames(table$Gene)
#
R <- table$r %>%
      matrix(nrow = 1) %>%
      magrittr::set_colnames(table$Gene)
#
samps <- BTF::collect_samples(samps_dir, last = 100) %>% BTF::align_samples()
#
gene_map <- readr::read_csv("gene_map.csv")
V <- apply(samps$V, c(2,3), mean) %>% magrittr::set_rownames(gene_map[["gene"]])
#
gns_to_use <- intersect(row.names(V), colnames(S))
V <- V[gns_to_use, ]
#
Ybar <- Ybar[, gns_to_use, drop = F] %>% magrittr::set_rownames("diff_score")
S <- S[, gns_to_use, drop = F] %>% magrittr::set_rownames("diff_score")
R <- R[, gns_to_use, drop = F] %>% magrittr::set_rownames("diff_score")
#
p_vals <- BTF::fabp_lin_reg(Y = Ybar, S = S, R = R, U = matrix(1, nrow = 1, ncol = 1), V = V, snr = 2)
#
wes_pal <- wesanderson::wes_palette("Zissou1", 5)
#
rank_p <- data.frame(p = p_vals$p, fdr = p_vals$fdr_p, rank = rank(p_vals$p), type = "UMPU", row = p_vals$row)
rank_fabp <- data.frame(p = p_vals$fabp, fdr = p_vals$fdr_fabp, rank = rank(p_vals$fabp), type = "FAB", row = p_vals$row)
rank_os <- data.frame(p = p_vals$p/2, fdr = p.adjust(p_vals$p/2, method = "BH"), rank = rank(p_vals$p), type = "OS", row = p_vals$row)
#
p1 <- rbind(rank_p, rank_fabp, rank_os) %>%
  dplyr::filter(fdr < 0.12) %>%
  ggplot2::ggplot() +
  geom_line(aes(x = rank, y = fdr, group = type, colour = type)) +
  geom_hline(yintercept = 0.1, size = 0.5, colour = "#a50f15", alpha = 0.5) +
  theme_light() +
  guides(colour = F) +
  scale_colour_manual(values = c("#969696", "#9e9ac8", "#000000"), name = "Test") +
  labs(x = "Rank", y = "p")
#
pos_df <- dplyr::filter(p_vals %>%
                          dplyr::mutate(x = rank(observed, ties.method = "first")), fdr_fabp < 0.1, observed > 0) %>%
          dplyr::arrange(observed)
neg_df <- dplyr::filter(p_vals %>%
                          dplyr::mutate(x = rank(observed, ties.method = "first")), fdr_fabp < 0.1, observed < 0) %>%
          dplyr::arrange(observed)
pl1 <- p_vals %>%
  dplyr::mutate(x = rank(observed, ties.method = "first")) %>%
  ggplot2::ggplot() +
  geom_bar(aes(x = x, y = observed, fill = fdr_fabp < 0.1), stat = "identity", colour = "white") +
  geom_text(data = pos_df, x = 100, y = seq(0.5, 2, length.out = nrow(pos_df)), aes(label = column), size = 3, hjust = 1.1, fontface = "bold") +
  geom_text(data = neg_df, x = 25, y = seq(-2, -0.5, length.out = nrow(neg_df)), aes(label = column), size = 3, hjust = -0.1, fontface = "bold") +
  geom_segment(data = pos_df, aes(x = x, y = observed), xend = 100, yend = seq(0.5, 2, length.out = nrow(pos_df)), colour = wes_pal[5], alpha = 0.5) +
  geom_segment(data = neg_df, aes(x = x, y = observed), xend = 25, yend = seq(-2, -0.5, length.out = nrow(neg_df)), colour = wes_pal[5], alpha = 0.5) +
  theme_light() +
  scale_fill_manual(values = c("#969696", wes_pal[3])) +
  ylim(c(-2, 2)) +
  guides(fill = F) +
  labs(x = "Genes ranked by differential CS", y = "Differential CS")
pos_df <- dplyr::filter(p_vals %>%
                          dplyr::mutate(x = rank(observed, ties.method = "first")), fdr_p < 0.1, observed > 0) %>%
  dplyr::arrange(observed)
neg_df <- dplyr::filter(p_vals %>%
                          dplyr::mutate(x = rank(observed, ties.method = "first")), fdr_p < 0.1, observed < 0) %>%
  dplyr::arrange(observed)
pl2 <- p_vals %>%
  dplyr::mutate(x = rank(observed, ties.method = "first")) %>%
  ggplot2::ggplot() +
  geom_bar(aes(x = x, y = observed, fill = fdr_p < 0.1), stat = "identity", colour = "white") +
  geom_text(data = pos_df, x = 100, y = seq(0.5, 2, length.out = nrow(pos_df)), aes(label = column), size = 3, hjust = 1.1, fontface = "bold") +
  geom_text(data = neg_df, x = 25, y = seq(-2, -0.5, length.out = nrow(neg_df)), aes(label = column), size = 3, hjust = -0.1, fontface = "bold") +
  geom_segment(data = pos_df, aes(x = x, y = observed), xend = 100, yend = seq(0.5, 2, length.out = nrow(pos_df)), colour = wes_pal[5], alpha = 0.5) +
  geom_segment(data = neg_df, aes(x = x, y = observed), xend = 25, yend = seq(-2, -0.5, length.out = nrow(neg_df)), colour = wes_pal[5], alpha = 0.5) +
  theme_light() +
  scale_fill_manual(values = c("#969696", wes_pal[2])) +
  ylim(c(-2, 2)) +
  guides(fill = F) +
  labs(x = "Genes ranked by differential CS", y = "Differential CS")
#
cowplot::plot_grid(pl1, pl2, nrow = 2, align = "hv")
ggplot2::ggsave("../../figs/wang_focused_example.pdf", height = 10, width = 10)

