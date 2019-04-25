#'
#'
#'
collate_tensor <- function(matrix_list, fill = NA, names = paste0("matrix", 0:(length(matrix_list)-1))){

  require(plyr)
  require(magrittr)

  all_colnames <- plyr::llply(matrix_list, colnames) %>%
    unlist() %>%
    unique()

  all_rownames <- plyr::llply(matrix_list, row.names) %>%
    unlist() %>%
    unique()

  plyr::llply(matrix_list, function(m){

    matrix_template <- matrix(fill, nrow = length(all_rownames), ncol = length(all_colnames)) %>%
                        magrittr::set_rownames(all_rownames) %>%
                        magrittr::set_colnames(all_colnames)
    matrix_template[row.names(m), colnames(m)] <- m
    matrix_template

  }) %>%
    magrittr::set_names(names)

}
