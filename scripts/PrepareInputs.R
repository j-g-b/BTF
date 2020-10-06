require(magrittr)
require(tidyverse)
#
cl_info <- readr::read_csv("data/sample_info.csv")
#
BTF::prepare_ge_rds("data/CCLE_expression.csv", "data/rds_obj/depmap_gene_expression.rds")
BTF::prepare_crispr_rds("data/Achilles_gene_effect.csv", "data/rds_obj/depmap_crispr.rds")
BTF::prepare_rnai_rds("data/D2_combined_gene_dep_scores.csv", "data/rds_obj/depmap_rnai.rds", cl_info)
BTF::prepare_mut_rds("data/CCLE_mutations.csv", "data/rds_obj/depmap_mut.rds")
BTF::prepare_crispr_rds("data/Achilles_gene_dependency.csv", "data/rds_obj/depmap_dprob.rds")
