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
p1 <- list()
p2 <- list()
p3 <- list()
#
wes_pal <- wesanderson::wes_palette("Zissou1", 5)
#
contrasts <- matrix(c("ML", "LP", "ML", "MS", "MS", "LP"), nrow = 3)
#
for(k in 1:3){
  Y <- (rowMeans(y$E[, ct == contrasts[k, 1]]) - rowMeans(y$E[, ct == contrasts[k, 2]])) %>%
    matrix(ncol = 1) %>%
    magrittr::set_rownames(y$genes$SYMBOL) %>%
    magrittr::set_colnames("diff_score")
  #
  S <- (apply(y$E[, ct == contrasts[k, 1]] - y$E[, ct == contrasts[k, 2]], 1, sd)) %>%
    matrix(ncol = 1) %>%
    magrittr::set_rownames(y$genes$SYMBOL) %>%
    magrittr::set_colnames("diff_score")
  #
  R <- rep(3, nrow(Y)) %>% 
    matrix(ncol = 1) %>%
    magrittr::set_rownames(y$genes$SYMBOL) %>%
    magrittr::set_colnames("diff_score")
  #
  gene_map <- readr::read_csv("gene_map.csv")
  V <- apply(samps$V, c(2,3), mean) %>% magrittr::set_rownames(gene_map[["gene"]])
  gns_to_use <- intersect(row.names(V), row.names(S))
  V <- V[gns_to_use, ]
  #
  p_vals <- BTF::fabp_lin_reg(Y = t(Y[gns_to_use, , drop=F]), 
                              S = t(S[gns_to_use, , drop=F]), 
                              R = t(R[gns_to_use, , drop=F]), V = V, U = matrix(1, nrow = 1, ncol = 1), snr = 1)
  #
  rank_p <- data.frame(p = p_vals$p, fdr = p_vals$fdr_p, rank = rank(p_vals$p), type = "UMPU", observed = p_vals$statistic)
  rank_fabp <- data.frame(p = p_vals$fabp, fdr = p_vals$fdr_fabp, rank = rank(p_vals$fabp), type = "FAB", observed = p_vals$statistic)
  rank_os <- data.frame(p = p_vals$p/2, fdr = p.adjust(p_vals$p/2, method = "BH"), rank = rank(p_vals$p), type = "OS", observed = p_vals$statistic)
  if(k==1){
    p1[[paste0(k, "")]] <- rbind(rank_p, rank_fabp) %>%
      dplyr::filter(fdr < 0.12) %>%
      ggplot2::ggplot() +
      geom_line(aes(x = rank, y = fdr, group = type, colour = type)) +
      geom_hline(yintercept = 0.1, size = 0.5, colour = "#a50f15", alpha = 0.5) +
      theme_light() +
      scale_colour_manual(values = c(wes_pal[2], wes_pal[3]), name = "Test") +
      labs(x = "Rank", y = "FDR") +
      theme(legend.position = c(0.2, 0.5))
  } else {
    p1[[paste0(k, "")]] <- rbind(rank_p, rank_fabp) %>%
      dplyr::filter(fdr < 0.12) %>%
      ggplot2::ggplot() +
      geom_line(aes(x = rank, y = fdr, group = type, colour = type)) +
      geom_hline(yintercept = 0.1, size = 0.5, colour = "#a50f15", alpha = 0.5) +
      theme_light() +
      scale_colour_manual(values = c(wes_pal[2], wes_pal[3]), name = "Test") +
      labs(x = "Rank", y = "FDR") +
      guides(colour = F)
  }
  #
  p2[[paste0(k, "")]] <- p_vals %>%
    ggplot2::ggplot() +
    geom_hex(aes(y = predicted, x = observed), bins = 100) +
    theme_light() + 
    scale_fill_viridis_c(option = "magma") +
    guides(fill = F) +
    labs(x = "Observed", y = "Predicted")
  #
  p3[[paste0(k, "")]] <- rbind(rank_p, rank_fabp) %>%
    ggplot2::ggplot() +
    geom_point(aes(y = -log10(fdr), x = observed, colour = type), alpha = 0.5) +
    geom_line(data = rank_os %>% dplyr::select(-type), 
              aes(y = -log10(fdr), x = observed)) +
    scale_colour_manual(values = c(wes_pal[2], wes_pal[3]), name = "Test") +
    theme_light() +
    facet_wrap(~type) +
    labs(x = "Observed differential expression", y = "-log10(FDR)") +
    guides(colour = F)
  #
}
#
cowplot::plot_grid(cowplot::plot_grid(p1$`1`, p2$`1`, ncol = 1),
                   cowplot::plot_grid(p1$`2`, p2$`2`, ncol = 1),
                   cowplot::plot_grid(p1$`3`, p2$`3`, ncol = 1),
                   ncol = 3)
ggplot2::ggsave("../../figs/limma_example.pdf", width = 12, height = 8)      

