#'
#'
#'
collect_samples <- function(save_dir, last = NULL){
  #
  require(plyr)
  require(magrittr)
  #
  sample_list <- list()
  #
  dir_files <- dir(save_dir) %>%
    magrittr::extract(., order(as.numeric(gsub(".*_", "", gsub("[.]csv", "", .)))))
  if(!is.null(last)){
    file_nums <- as.numeric(gsub(".*_", "", gsub("[.]csv", "", dir_files)))
    max_iter <- max(file_nums)
    dir_files <- dir_files[file_nums > max_iter - last]
  }
  #
  if(any(grepl("^U_.*[.]csv", dir_files))){
    sample_list[["U"]] <- dir_files %>%
      magrittr::extract(., grepl("^U_.*[.]csv", .)) %>%
      paste0(save_dir, "/", .) %>%
      plyr::aaply(1, function(fn){

        read.csv(fn, header = F) %>%
          as.matrix()

      })
  }
  if(any(grepl("^V_.*[.]csv", dir_files))){
    sample_list[["V"]] <- dir_files %>%
      magrittr::extract(., grepl("^V_.*[.]csv", .)) %>%
      paste0(save_dir, "/", .) %>%
      plyr::aaply(1, function(fn){

        read.csv(fn, header = F) %>%
          as.matrix()

      })
  }
  if(any(grepl("^R.*_.*[.]csv", dir_files))){
    sample_list[["R"]] <- dir_files %>%
      magrittr::extract(., grepl("^R.*_.*[.]csv", .)) %>%
      data.frame() %>%
      magrittr::set_colnames("fname") %>%
      dplyr::mutate(rel_num = as.numeric(gsub("_.*", "", gsub("R", "", fname)))) %>%
      plyr::daply(.(rel_num), function(d){
        plyr::aaply(as.character(d[["fname"]]), 1, function(fn){
          read.csv(paste0(save_dir, "/", fn), header = F) %>% as.matrix()
        })
      }) %>%
      aperm(c(2, 1, 3, 4))
  }
  if(any(grepl("^SS_.*[.]csv", dir_files))){
    sample_list[["SS"]] <- dir_files %>%
      magrittr::extract(., grepl("^SS_.*[.]csv", .)) %>%
      paste0(save_dir, "/", .) %>%
      plyr::aaply(1, function(fn){

        read.csv(fn, header = F) %>%
          as.matrix() %>%
          c()

      })
  }
  if(any(grepl("^LambdaU_.*[.]csv", dir_files))){
    sample_list[["LambdaU"]] <- dir_files %>%
      magrittr::extract(., grepl("^LambdaU_.*[.]csv", .)) %>%
      paste0(save_dir, "/", .) %>%
      plyr::aaply(1, function(fn){

        read.csv(fn, header = F) %>%
          as.matrix()

      })
  }
  if(any(grepl("^LambdaV_.*[.]csv", dir_files))){
    sample_list[["LambdaV"]] <- dir_files %>%
      magrittr::extract(., grepl("^V_.*[.]csv", .)) %>%
      paste0(save_dir, "/", .) %>%
      plyr::aaply(1, function(fn){

        read.csv(fn, header = F) %>%
          as.matrix()

      })
  }
  if(any(grepl("^MuU_.*[.]csv", dir_files))){
    sample_list[["MuU"]] <- dir_files %>%
      magrittr::extract(., grepl("^MuU_.*[.]csv", .)) %>%
      paste0(save_dir, "/", .) %>%
      plyr::aaply(1, function(fn){

        read.csv(fn, header = F) %>%
          as.matrix() %>%
          c()

      })
  }
  if(any(grepl("^MuV_.*[.]csv", dir_files))){
    sample_list[["MuV"]] <- dir_files %>%
      magrittr::extract(., grepl("^MuV_.*[.]csv", .)) %>%
      paste0(save_dir, "/", .) %>%
      plyr::aaply(1, function(fn){

        read.csv(fn, header = F) %>%
          as.matrix() %>%
          c()

      })
  }
  if(any(grepl("^XiR_.*[.]csv", dir_files))){
    sample_list[["XiR"]] <- dir_files %>%
      magrittr::extract(., grepl("^XiR_.*[.]csv", .)) %>%
      paste0(save_dir, "/", .) %>%
      plyr::aaply(1, function(fn){

        read.csv(fn, header = F) %>%
          as.matrix() %>%
          c()

      })
  }
  if(any(grepl("^PsiR_.*[.]csv", dir_files))){
    sample_list[["PsiR"]] <- dir_files %>%
      magrittr::extract(., grepl("^PsiR_.*[.]csv", .)) %>%
      paste0(save_dir, "/", .) %>%
      plyr::aaply(1, function(fn){

        read.csv(fn, header = F) %>%
          as.matrix()

      })
  }
  return(sample_list)
}
