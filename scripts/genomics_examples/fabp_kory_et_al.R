#
samps_dir <- "../../alt_samps"
#
library(tidyverse)
library(magrittr)
#
K562_table <- readr::read_csv("K562 SHMT1-null.csv") %>%
                dplyr::filter(!grepl("INTERGENIC", sgRNA)) %>%
                dplyr::mutate(gene = gsub("_.*", "", gsub("sg", "", sgRNA)),
                              replicate = as.integer(gsub(".*_", "", sgRNA)),
                              full_lfc = log2(((`full media` + 1)/sum(`full media` + 1, na.rm = T)) / ((`initial` + 1)/sum(`initial` + 1, na.rm = T))),
                              minus_lfc = log2(((`minus serine` + 1)/sum(`minus serine` + 1, na.rm = T)) / ((`initial` + 1)/sum(`initial` + 1, na.rm = T)))) %>%
                dplyr::select(-sgRNA) %>%
                dplyr::group_by(gene) %>%
                dplyr::summarise(r = sum(!is.na(minus_lfc - full_lfc)),
                                 s = sd(minus_lfc - full_lfc, na.rm = T),
                                 m = mean(minus_lfc - full_lfc, na.rm = T)) %>%
                dplyr::filter(r > 1)
#
Jurkat_table <- readr::read_csv("Jurkat SHMT1-null.csv") %>%
                  dplyr::filter(!grepl("INTERGENIC", sgRNA)) %>%
                    dplyr::mutate(gene = gsub("_.*", "", gsub("sg", "", sgRNA)),
                                  replicate = as.integer(gsub(".*_", "", sgRNA)),
                                  full_lfc = log2(((`full media` + 1)/sum(`full media` + 1, na.rm = T)) / ((`initial` + 1)/sum(`initial` + 1, na.rm = T))),
                                  minus_lfc = log2(((`minus serine` + 1)/sum(`minus serine` + 1, na.rm = T)) / ((`initial` + 1)/sum(`initial` + 1, na.rm = T)))) %>%
                  dplyr::select(-sgRNA) %>%
                  dplyr::group_by(gene) %>%
                  dplyr::summarise(r = sum(!is.na(minus_lfc - full_lfc)),
                                   s = sd(minus_lfc - full_lfc, na.rm = T),
                                   m = mean(minus_lfc - full_lfc, na.rm = T)) %>%
                  dplyr::filter(r > 1)
#
Ybar <- rbind(data.frame(cell_line = "ACH-000995", 
                         gene = Jurkat_table[["gene"]],
                         value = Jurkat_table[["m"]]),
              data.frame(cell_line = "ACH-000551", 
                         gene = K562_table[["gene"]],
                         value = K562_table[["m"]])) %>%
        reshape2::acast(cell_line ~ gene, value.var = "value")
cols_to_use <- apply(Ybar, 2, function(x){!any(is.na(x))})
Ybar <- Ybar[, cols_to_use]
#
S <- rbind(data.frame(cell_line = "ACH-000995", 
                         gene = Jurkat_table[["gene"]],
                         value = Jurkat_table[["s"]]),
              data.frame(cell_line = "ACH-000551", 
                         gene = K562_table[["gene"]],
                         value = K562_table[["s"]])) %>%
     reshape2::acast(cell_line ~ gene, value.var = "value")
S <- S[, cols_to_use]
#
R <- rbind(data.frame(cell_line = "ACH-000995", 
                      gene = Jurkat_table[["gene"]],
                      value = Jurkat_table[["r"]]),
           data.frame(cell_line = "ACH-000551", 
                      gene = K562_table[["gene"]],
                      value = K562_table[["r"]])) %>%
     reshape2::acast(cell_line ~ gene, value.var = "value")
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
p_vals <- BTF::fabp_lin_reg(Ybar, S, R, U = U, V = V)
#
wes_pal <- wesanderson::wes_palette("Zissou1", 5)
#
pl1 <- p_vals %>%
  ggplot2::ggplot() +
  geom_point(aes(x = fdr_p, y = fdr_fabp, 
                 colour = ifelse(fdr_fabp <= 0.1 & fdr_p > 0.1, wes_pal[3], ifelse(fdr_fabp > 0.1 & fdr_p <= 0.1, wes_pal[2], "#969696"))), 
             alpha = 0.75) +
  geom_abline(slope = 1, colour = wes_pal[5]) +
  theme_light() +
  scale_colour_manual(values = c(wes_pal[2], "#969696", wes_pal[3]), name = "", labels = c("Only UMPU", "Neither", "Only FAB")) +
  theme(legend.position = c(0.8, 0.2)) +
  labs(x = "FDR UMPU", y = "FDR FAB")
#
pl2 <- p_vals %>%
  ggplot2::ggplot() +
  geom_hex(aes(y = predicted, x = observed), bins = 100) +
  theme_light() + 
  scale_fill_viridis_c(option = "magma") +
  guides(fill = F) +
  labs(x = "Observed", y = "Predicted")
#
prop_common_discov_1 <- rep(0, 100 - 4)
for(n in 5:100){
  k <- p_vals %>% 
        dplyr::filter(row == "ACH-000551") %>%
        dplyr::arrange(fdr_p, column) %>%
        dplyr::mutate(rank = 1:n()) %>%
        dplyr::filter(rank <= n) %>%
        magrittr::extract2("column")
  j <- p_vals %>% 
    dplyr::filter(row == "ACH-000551") %>%
    dplyr::arrange(fdr_fabp, column) %>%
    dplyr::mutate(rank = 1:n()) %>%
    dplyr::filter(rank <= n) %>%
    magrittr::extract2("column")
  common_discov <- length(intersect(j, k)) / n
  prop_common_discov_1[n - 4] <- common_discov
}
#
prop_common_discov_2 <- rep(0, 100 - 4)
for(n in 5:100){
  k <- p_vals %>% 
    dplyr::filter(row == "ACH-000995") %>%
    dplyr::arrange(fdr_p, column) %>%
    dplyr::mutate(rank = 1:n()) %>%
    dplyr::filter(rank <= n) %>%
    magrittr::extract2("column")
  j <- p_vals %>% 
    dplyr::filter(row == "ACH-000995") %>%
    dplyr::arrange(fdr_fabp, column) %>%
    dplyr::mutate(rank = 1:n()) %>%
    dplyr::filter(rank <= n) %>%
    magrittr::extract2("column")
  common_discov <- length(intersect(j, k)) / n
  prop_common_discov_2[n - 4] <- common_discov
}
#
pl3 <- data.frame(prop = c(prop_common_discov_1, prop_common_discov_2),
           cell_line = c(rep("ACH-000551", length(prop_common_discov_1)),
                         rep("ACH-000995", length(prop_common_discov_2))),
           total_discoveries = c(5:100, 5:100)) %>%
  ggplot2::ggplot() +
  geom_area(aes(x = total_discoveries, y = prop), alpha = 0.6) +
  facet_wrap(~cell_line, nrow = 2) +
  theme_light() +
  ylim(c(0, 1)) +
  labs(y = "Proportion shared discoveries", x = "Total discoveries")
#
pl4 <- cowplot::plot_grid(pl2, pl1)
cowplot::plot_grid(pl4, pl3, rel_widths = c(2, 1))
#
ggplot2::ggsave("../../figs/kory_example.pdf", height = 6, width = 16)
