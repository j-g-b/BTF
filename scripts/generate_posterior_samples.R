library(tidyverse)
library(magrittr)
#
set.seed(59)
#
ge <- readRDS("data/rds_obj/depmap_gene_expression.rds")
cr <- readRDS("data/rds_obj/depmap_crispr.rds")
rnai <- readRDS("data/rds_obj/depmap_rnai.rds")
mut <- readRDS("data/rds_obj/depmap_mut.rds")
cl_info <- readr::read_csv("data/sample_info.csv")
#
ge %<>% magrittr::set_rownames(., gsub(" [(].*", "", row.names(.)))
cr %<>% magrittr::set_rownames(., gsub(" [(].*", "", row.names(.)))
rnai %<>% magrittr::set_rownames(., gsub(" [(].*", "", row.names(.)))
#
ge <- ge[!duplicated(row.names(ge)), !duplicated(colnames(ge))]
cr <- cr[!duplicated(row.names(cr)), !duplicated(colnames(cr))]
rnai <- rnai[!duplicated(row.names(rnai)), !duplicated(colnames(rnai))]
mut <- mut[!duplicated(row.names(mut)), !duplicated(colnames(mut))]
#
dim1 <- max(ncol(ge), ncol(cr), ncol(rnai), ncol(mut))
#
top_2m_exp_genes <- apply(ge, 1, function(x){mean(x^2)}) %>%
                    rank(1/.) %>%
                    data.frame(gene = names(.), rank = .) %>%
                    dplyr::filter(rank < 1001)
top_var_exp_genes <- apply(ge, 1, function(x){var(x)}) %>%
                      rank(1/.) %>%
                      data.frame(gene = names(.), rank = .) %>%
                      dplyr::filter(rank < 1001)
census_genes <- readr::read_csv("data/cancer_gene_census.csv")[["Gene Symbol"]] %>%
                  c("TLCD1") %>%
                  union(readr::read_csv("data/wang_focused_genes.csv")[["gene"]]) %>%
                  union(readr::read_csv("data/kory_genes.csv")[["gene"]]) %>%
                  union(top_2m_exp_genes[["gene"]]) %>%
                  union(top_var_exp_genes[["gene"]]) %>%
                  intersect(row.names(ge))
cls <- cl_info %>%
        dplyr::filter(!(`disease` %in% c("Fibroblast", "fibroblast", "Immortalized"))) %>%
        magrittr::extract2("DepMap_ID") %>%
        intersect(colnames(ge))
#
N <- cls
M <- length(census_genes)
d1 <- 16
d2 <- 64
S <- 5000
#
TensorList <- BTF::collate_tensor(list(t(ge), t(cr), t(rnai), t(mut))) %>%
              BTF::subset_tensor(cls, census_genes)
#
K <- length(TensorList)
X <- BTF::run_btf(n_samples = S, tensor_list = TensorList,
                  d1 = d1, d2 = d2, save_dir = "alt_samps",
                  burn = 1000, thin = 40, matrix_type = c(2, 0, 0, 1))
