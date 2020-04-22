#
samps_dir <- "../../alt_samps"
#
library(tidyverse)
library(magrittr)
#
treatment_info <- readr::read_csv("primary_replicate_treatment_info.csv")
table <- readr::read_csv("primary_logfold_change.csv") %>%
  dplyr::rename(cl = X1) %>%
  reshape2::melt(id.vars = "cl", variable.name = "column_name") %>%
  dplyr::left_join(treatment_info) %>%
  dplyr::group_by(broad_id, cl) %>%
  dplyr::filter(length(unique(dose)) == 1) %>%
  dplyr::filter(perturbation_type == "experimental_treatment") %>%
  dplyr::ungroup()
#
Ybar <- table %>%
  reshape2::acast(cl ~ broad_id, 
                  value.var = "value", 
                  fun.aggregate = function(x){if(length(x) == 1){NA}else{mean(x, na.rm = T)}})
#
S <- table %>%
  reshape2::acast(cl ~ broad_id, 
                  value.var = "value",
                  fun.aggregate = function(x){if(length(x) == 1){NA}else{sd(x, na.rm = T)}})
#
R <- table %>%
  reshape2::acast(cl ~ broad_id, 
                  value.var = "value",
                  fun.aggregate = function(x){sum(!is.na(x))})
#
samps <- BTF::collect_samples(samps_dir, last = 100) %>% BTF::align_samples()
#
cl_map <- readr::read_csv("cl_map.csv")
U <- apply(samps$U, c(2,3), mean) %>% magrittr::set_rownames(cl_map[["depmap_id"]])
#
all_cpd_df <- data.frame()
for(drug_name in colnames(Ybar)){
  cls_to_use <- intersect(row.names(U), row.names(S)[!is.na(S[, drug_name])])
  Ud <- U[cls_to_use, ]
  #
  Ybard <- Ybar[cls_to_use, drug_name, drop = F]
  Sd <- S[cls_to_use, drug_name, drop = F]
  Rd <- R[cls_to_use, drug_name, drop = F]
  #
  p_vals <- BTF::fabp_lin_reg(Ybard, Sd, Rd, Ud, V = matrix(1, nrow = 1, ncol = 1), snr = 2)
  all_cpd_df <- rbind(all_cpd_df,
                      data.frame(n_discov_fab = c(sum(p_vals$fdr_fabp < 0.1), sum(p_vals$fdr_fabp < 0.05), sum(p_vals$fdr_fabp < 0.01), sum(p_vals$fdr_fabp < 0.001)),
                                 n_discov_p = c(sum(p_vals$fdr_p < 0.1), sum(p_vals$fdr_p < 0.05), sum(p_vals$fdr_p < 0.01), sum(p_vals$fdr_p < 0.001)),
                                 n_discov_os = c(sum(p.adjust(p_vals$p/2, method = "BH") < 0.1), sum(p.adjust(p_vals$p/2, method = "BH") < 0.05), sum(p.adjust(p_vals$p/2, method = "BH") < 0.01), sum(p.adjust(p_vals$p/2, method = "BH") < 0.001)),
                                 fdr = c(0.1, 0.05, 0.01, 0.001),
                                 n_cls = nrow(p_vals),
                                 r2 = p_vals$r2[1],
                                 mv = p_vals$model_error[1],
                                 ev = p_vals$error_var[1]))
  #
  cat(paste0("\r", drug_name))
}
#
wes_pal <- wesanderson::wes_palette("Zissou1", 5)
#
pdf <- all_cpd_df %>%
  dplyr::mutate(n_discov_fab = ifelse(n_discov_os > 0, (n_discov_fab) / (n_discov_os), 1),
                n_discov_p = ifelse(n_discov_os > 0, (n_discov_p) / (n_discov_os), 1)) %>%
  dplyr::group_by(n_discov_fab, n_discov_p, fdr) %>%
  dplyr::summarise(count = log2(n())) %>%
  dplyr::ungroup() 

pl1df <- all_cpd_df %>% 
        dplyr::filter(fdr == 0.1) %>%
        dplyr::mutate(n_discov_diff = n_discov_fab - n_discov_p) %>%
        dplyr::mutate(x = rank(n_discov_diff, ties.method = "first"))
pl1 <- pl1df %>% 
        ggplot2::ggplot() + 
        geom_abline(slope = 1, colour = wes_pal[5]) +
        geom_point(aes(x = n_discov_p, y = n_discov_fab, 
                       colour = ifelse(n_discov_diff > 0, "a", ifelse(n_discov_diff == 0, "b", "c"))),
                   size = 1, shape = 1) + 
        theme_light() + 
        scale_colour_manual(values = c(wes_pal[3], "#969696", wes_pal[2]), name = "", labels = c("More discoveries", 
                                                                                                 "Equal discoveries",
                                                                                                 "Fewer discoveries")) +
        labs(x = "Number of Discoveries UMPU", y = "Number of Discoveries FAB") +
        theme(legend.position = c(0.2, 0.8))
pl1
ggplot2::ggsave("../../figs/repurposing_example2.pdf", height = 5, width = 5)