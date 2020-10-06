#
samps_dir <- "../../alt_samps"
#
library(tidyverse)
library(magrittr)
#
samps <- BTF::collect_samples(samps_dir, last = 100) %>% BTF::align_samples()
#
n <- 12
r <- n/2
spike_props <- c(0.5, 0.25, 0.1, 0.05)
#
# Load Sanger data
sanger_dep <- readRDS("../../../data/rds_obj/sanger_crispr.rds") %>% magrittr::set_rownames(., gsub(" [(].*", "", row.names(.)))
#
gene_map <- readr::read_csv("gene_map.csv")
V <- apply(samps$V, c(2,3), mean) %>% magrittr::set_rownames(gene_map[["gene"]])
gns_to_use <- intersect(row.names(V), row.names(sanger_dep))
V <- V[gns_to_use, ]
sanger_dep <- sanger_dep[gns_to_use, ]
#
plist <- list()
#
for(sp in 1:length(spike_props)){
  #
  spike_prop <- spike_props[sp]
  # Spike-in signal from N(0, something); keep track of where the signals are spiked-in
  spike_in <- rbernoulli(nrow(V), p = spike_prop)
  print(sum(spike_in) / nrow(V))
  spike_effects <- rnorm(sum(spike_in), sd = 0.25)
  #
  p_vals <- list()
  #
  for(t in 1:1000){
    #
    print(t)
    # Randomly sample n cell lines; randomly assign them to one of two groups
    r_cls <- sample(colnames(sanger_dep), n)
    r_cls_a <- sample(r_cls, r)
    r_cls_b <- setdiff(r_cls, r_cls_a)
    #
    Y <- rowMeans(apply(sanger_dep[, r_cls_a], 2, function(c){c[spike_in] <- c[spike_in] + spike_effects/2; c})) %>%
      matrix(nrow = 1) %>%
      magrittr::set_colnames(row.names(V)) %>%
      magrittr::set_rownames("score")
    #
    Y1 <- rowMeans(apply(sanger_dep[, r_cls_b], 2, function(c){c[spike_in] <- c[spike_in] - spike_effects/2; c})) %>%
      matrix(nrow = 1) %>%
      magrittr::set_colnames(row.names(V)) %>%
      magrittr::set_rownames("score")
    #
    S <- apply(sanger_dep[, r_cls_a], 1, sd) %>%
      matrix(nrow = 1) %>%
      magrittr::set_colnames(row.names(V)) %>%
      magrittr::set_rownames("score")
    #
    S1 <- apply(sanger_dep[, r_cls_b], 1, sd) %>%
      matrix(nrow = 1) %>%
      magrittr::set_colnames(row.names(V)) %>%
      magrittr::set_rownames("score")
    #
    pval_df <- BTF::fabp_lin_reg(Y = Y, S = S, R = matrix(r, 1, nrow(V)),
                                 U = matrix(1, nrow = 1, ncol = 1), V = V,
                                 pool_sampling_var = F,
                                 Y1 = Y1, S1 = S1, R1 = matrix(r, 1, nrow(V)))
    p_vals[[paste0("", t)]] <- pval_df %>% dplyr::mutate(trial = t, spiked = spike_in)
  }
  #
  p_val_df <- data.table::rbindlist(p_vals)
  #
  wes_pal <- wesanderson::wes_palette("Zissou1", 5)
  #
  plist[[sp]] <- p_val_df %>%
    dplyr::group_by(trial) %>%
    dplyr::summarise(UMPU = sum(fdr_p < 0.1 & !spiked) / max(1, sum(fdr_p < 0.1)),
                     FAB = sum(fdr_fabp < 0.1 & !spiked) / max(1, sum(fdr_fabp < 0.1))) %>%
    dplyr::ungroup() %>%
    reshape2::melt(id.vars = "trial", variable.name = "type") %>%
    dplyr::group_by(type) %>%
    dplyr::summarise(mc_mean = mean(value),
                     mc_sd = sd(value) / sqrt(n())) %>%
    ggplot2::ggplot() +
    geom_bar(aes(x = type, y = mc_mean, fill = type), stat = "identity") +
    geom_errorbar(aes(x = type, ymin = mc_mean - mc_sd, ymax = mc_mean + mc_sd), width = 0.2) +
    geom_hline(yintercept = 0.1, colour = wes_pal[5]) +
    geom_hline(yintercept = 0.1*(1 - spike_prop), colour = wes_pal[5], linetype = 2) +
    scale_fill_manual(values = c(wes_pal[2], wes_pal[4]), name = "Type") +
    theme_light() +
    guides(fill = F) +
    labs(x = "Test", y = "Empirical FDR", title = paste0(100*spike_prop, "% non-null"))
}
#
cowplot::plot_grid(plotlist = plist)
ggplot2::ggsave("../../figs/spike_in_simulation_study_plot.pdf", width = 7, height = 7)
