#
samps_dir <- "../../alt_samps"
#
library(tidyverse)
library(magrittr)
#
table <- readr::read_tsv("pool.normalized.counts.txt") %>%
  dplyr::filter(!grepl("INTERGENIC|CTRL", sgRNA)) %>%
  reshape2::melt(id.vars = "sgRNA", variable.name = "cell_cond") %>%
  dplyr::mutate(cl = gsub("-final|-initial", "", cell_cond)) %>%
  dplyr::mutate(fi = ifelse(grepl("final", cell_cond), "final", "initial")) %>%
  dplyr::select(-cell_cond) %>%
  reshape2::dcast(sgRNA + cl ~ fi) %>%
  dplyr::filter(!is.na(final), !is.na(initial)) %>%
  dplyr::mutate(gene = gsub("_.*", "", gsub("sg", "", sgRNA)),
                replicate = as.integer(gsub(".*_", "", sgRNA))) %>%
  dplyr::group_by(cl) %>%
  dplyr::mutate(lfc = log2(((`final` + 1)/sum(`final` + 1, na.rm = T)) / ((`initial` + 1)/sum(`initial` + 1, na.rm = T)))) %>%
  dplyr::ungroup() %>%
  dplyr::select(-sgRNA) %>%
  dplyr::group_by(gene, cl) %>%
  dplyr::summarise(r = sum(!is.na(lfc)),
                   s = sd(lfc, na.rm = T),
                   m = mean(lfc, na.rm = T)) %>%
  dplyr::filter(r > 1)
#
Ybar <- table %>%
  dplyr::left_join(readr::read_csv("cl_conversion.csv")) %>%
  dplyr::select(-cl) %>%
  reshape2::acast(depmap_id ~ gene, value.var = "m")
cols_to_use <- apply(Ybar, 2, function(x){!any(is.na(x))})
Ybar <- Ybar[, cols_to_use]
#
S <- table %>%
  dplyr::left_join(readr::read_csv("cl_conversion.csv")) %>%
  dplyr::select(-cl) %>%
  reshape2::acast(depmap_id ~ gene, value.var = "s")
S <- S[, cols_to_use]
#
R <- table %>%
  dplyr::left_join(readr::read_csv("cl_conversion.csv")) %>%
  dplyr::select(-cl) %>%
  reshape2::acast(depmap_id ~ gene, value.var = "r")
R <- R[, cols_to_use]
#
samps <- BTF::collect_samples(samps_dir, last = 100) %>% BTF::align_samples()
#
cl_map <- readr::read_csv("cl_map.csv")
gene_map <- readr::read_csv("gene_map.csv")
V <- apply(samps$V, c(2,3), mean) %>% magrittr::set_rownames(gene_map[["gene"]])
U <- apply(samps$U, c(2,3), mean) %>% magrittr::set_rownames(cl_map[["depmap_id"]])
#
cls_to_use <- intersect(row.names(U), row.names(S))
gns_to_use <- intersect(row.names(V), colnames(S))
V <- V[gns_to_use, ]
U <- U[cls_to_use, ]
#
Ybar <- Ybar[cls_to_use, gns_to_use]
S <- S[cls_to_use, gns_to_use]
R <- R[cls_to_use, gns_to_use]
#
p_vals <- BTF::fabp_lin_reg(Y = Ybar, S = S, R = R, U = U, V = V, pool_sampling_var = F)
#
cl_conv <- readr::read_csv("cl_conversion.csv")
cls_with_cr <- cl_conv$depmap_id[cl_conv$has_crispr == 1]
cls_witho_cr <- cl_conv$depmap_id[cl_conv$has_crispr == 0]
#
wes_pal <- wesanderson::wes_palette("Zissou1", 5)
#
rank_p <- data.frame(p = p_vals$p, fdr = p_vals$fdr_p, rank = rank(p_vals$p), type = "UMPU", row = p_vals$row)
rank_fabp <- data.frame(p = p_vals$fabp, fdr = p_vals$fdr_fabp, rank = rank(p_vals$fabp), type = "FAB", row = p_vals$row)
#
p1 <- rbind(rank_p, rank_fabp) %>%
  dplyr::mutate(type = factor(type, levels = c("UMPU", "FAB"))) %>%
  dplyr::filter(fdr < 0.12) %>%
  ggplot2::ggplot() +
  geom_line(aes(x = rank, y = fdr, group = type, colour = type)) +
  geom_hline(yintercept = 0.1, size = 0.5, colour = wes_pal[5], alpha = 0.5) +
  theme_light() +
  guides(colour = F) +
  scale_colour_manual(values = c(wes_pal[1], wes_pal[4]), name = "Test") +
  labs(x = "Rank", y = "FDR")
#
p2 <- p_vals %>%
  dplyr::filter(row %in% cls_with_cr) %>%
  ggplot2::ggplot() +
  geom_hex(aes(y = predicted, x = observed), bins = 200) +
  theme_light() +
  scale_fill_viridis_c(option = "magma") +
  guides(fill = F) +
  labs(x = "Observed", y = "Predicted")
#
p3 <- p_vals %>%
  dplyr::filter(row %in% cls_witho_cr) %>%
  ggplot2::ggplot() +
  geom_hex(aes(y = predicted, x = observed), bins = 200) +
  theme_light() +
  scale_fill_viridis_c(option = "magma") +
  guides(fill = F) +
  labs(x = "Observed", y = "Predicted")
#
p4 <- rbind(rank_p, rank_fabp) %>%
  dplyr::mutate(type = factor(type, levels = c("UMPU", "FAB"))) %>%
  dplyr::filter(fdr < 0.1) %>%
  ggplot2::ggplot() +
  geom_bar(aes(x = type, fill = type)) +
  theme_light() +
  scale_fill_manual(values = c(wes_pal[1], wes_pal[4]), name = "Test") +
  guides(fill = F) +
  labs(x = "Test", y = "# of discoveries")
#
cowplot::plot_grid(p4, p1, p2, p3, align = "hv")
ggplot2::ggsave("../../figs/wang_example.pdf", width = 7, height = 6)
