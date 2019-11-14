#'
#'
#'
warm_start <- function(save_dir){
  #
  require(plyr)
  require(magrittr)
  #
  inits <- list()
  #
  dir_files <- dir(save_dir)
  #
  if(any(grepl("^U_.*[.]csv", dir_files))){
    inits[["U0"]] <- dir_files %>%
                      magrittr::extract(., grepl("^U_.*[.]csv", .)) %>%
                      magrittr::extract(., as.numeric(gsub(".*_", "", gsub("[.]csv", "", .))) == max(as.numeric(gsub(".*_", "", gsub("[.]csv", "", .))))) %>%
                      paste0(save_dir, "/", .) %>%
                      read.csv(header = F) %>%
                      as.matrix()
  }
  if(any(grepl("^V_.*[.]csv", dir_files))){
    inits[["V0"]] <- dir_files %>%
                      magrittr::extract(., grepl("^V_.*[.]csv", .)) %>%
                      magrittr::extract(., as.numeric(gsub(".*_", "", gsub("[.]csv", "", .))) == max(as.numeric(gsub(".*_", "", gsub("[.]csv", "", .))))) %>%
                      paste0(save_dir, "/", .) %>%
                      read.csv(header = F) %>%
                      as.matrix()
  }
  if(any(grepl("^R.*_.*[.]csv", dir_files))){
    inits[["R0"]] <- dir_files %>%
                      magrittr::extract(., grepl("^R.*_.*[.]csv", .)) %>%
                      magrittr::extract(., as.numeric(gsub(".*_", "", gsub("[.]csv", "", .))) == max(as.numeric(gsub(".*_", "", gsub("[.]csv", "", .))))) %>%
                      magrittr::extract(., order(as.numeric(gsub("_.*", "", gsub("R", "", .))))) %>%
                      paste0(save_dir, "/", .) %>%
                      plyr::alply(1, function(fn){
                        read.csv(fn, header = F) %>% as.matrix()
                      }) %>%
                      magrittr::set_names(., paste0("matrix", as.numeric(names(.)) - 1))
  }
  if(any(grepl("^SS_.*[.]csv", dir_files))){
    inits[["SS0"]] <- dir_files %>%
                      magrittr::extract(., grepl("^SS_.*[.]csv", .)) %>%
                      magrittr::extract(., as.numeric(gsub(".*_", "", gsub("[.]csv", "", .))) == max(as.numeric(gsub(".*_", "", gsub("[.]csv", "", .))))) %>%
                      paste0(save_dir, "/", .) %>%
                      read.csv(header = F) %>%
                      as.matrix() %>%
                      c()
  }
  if(any(grepl("^Mu_.*[.]csv", dir_files))){
    inits[["Mu0"]] <- dir_files %>%
      magrittr::extract(., grepl("^Mu_.*[.]csv", .)) %>%
      magrittr::extract(., as.numeric(gsub(".*_", "", gsub("[.]csv", "", .))) == max(as.numeric(gsub(".*_", "", gsub("[.]csv", "", .))))) %>%
      paste0(save_dir, "/", .) %>%
      read.csv(header = F) %>%
      as.matrix() %>%
      c()
  }
  inits[["IterStart"]] <- max(as.numeric(gsub(".*_", "", gsub("[.]csv", "", dir_files)))) + 1
  return(inits)
}
