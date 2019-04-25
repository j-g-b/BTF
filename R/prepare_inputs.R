#'
#'
#'
prepare_ge_rds <- function(input_path, output_path){

  res <- readr::read_csv(input_path) %>%
    magrittr::set_rownames(., magrittr::extract2(., "X1")) %>%
    dplyr::select(-X1) %>%
    as.matrix() %>%
    t() %>%
    saveRDS(file = output_path)

}

#'
#'
#'
prepare_crispr_rds <- function(input_path, output_path){

  readr::read_csv(input_path) %>%
    magrittr::set_rownames(., magrittr::extract2(., "X1")) %>%
    dplyr::select(-X1) %>%
    as.matrix() %>%
    t() %>%
    saveRDS(file = output_path)

}

#'
#'
#'
prepare_cn_rds <- function(input_path, output_path){

  readr::read_csv(input_path) %>%
    magrittr::set_rownames(., magrittr::extract2(., "X1")) %>%
    dplyr::select(-X1) %>%
    as.matrix() %>%
    t() %>%
    magrittr::raise_to_power(1.1, .) %>%
    saveRDS(file = output_path)

}

#'
#'
#'
prepare_rnai_rds <- function(input_path, output_path, cl_info){
  require(magrittr)
  cl_map <- cl_info %>%
              dplyr::select(DepMap_ID, CCLE_Name) %>%
              magrittr::set_rownames(., magrittr::extract2(., "CCLE_Name")) %>%
              dplyr::select(-CCLE_Name) %>%
              as.matrix()
  res <-readr::read_csv(input_path) %>%
    magrittr::set_rownames(., magrittr::extract2(., "X1")) %>%
    dplyr::select(-X1) %>%
    as.matrix()
  res <- res[, intersect(colnames(res), row.names(cl_map))]
  res %>%
    magrittr::set_colnames(., cl_map[intersect(colnames(.), row.names(cl_map)), 1]) %>%
    saveRDS(file = output_path)

}

#'
#'
#'
prepare_mut_rds <- function(input_path, output_path){
  require(magrittr)
  readr::read_csv(input_path) %>%
    dplyr::filter(Variant_annotation %in% c("damaging", "other non-conserving")) %>%
    dplyr::mutate(value = 1) %>%
    reshape2::acast(Hugo_Symbol~DepMap_ID, value.var = "value", fun.aggregate = function(x){length(unique(x))}) %>%
    saveRDS(file = output_path)

}

#'
#'
#'
prepare_rppa_rds <- function(input_path, output_path, cl_info){

}

#'
#'
#'
prepare_gdsc_auc_rds <- function(input_path, output_path){
  require(magrittr)
  readr::read_csv(input_path) %>%
    magrittr::set_rownames(., magrittr::extract2(., "X1")) %>%
    dplyr::select(-X1) %>%
    as.matrix() %>%
    saveRDS(file = output_path)
}

#'
#'
#'
prepare_gdsc_ic50_rds <- function(input_path, output_path){
  require(magrittr)
  readr::read_csv(input_path) %>%
    magrittr::set_rownames(., magrittr::extract2(., "X1")) %>%
    dplyr::select(-X1) %>%
    as.matrix() %>%
    saveRDS(file = output_path)
}
