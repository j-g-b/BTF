#'
#'
#'
evaluate_val_r2 <- function(validation_set, samples, tensor_list, type = "global"){
  #
  require(plyr)
  require(magrittr)
  #
  n_samps <- samples[["U"]] %>% nrow()
  #
  if(type == "global"){
    plyr::alply(1:length(validation_set), 1, function(v){
      val_df <- validation_set[[v]]
      plyr::aaply(1:n_samps, 1, function(s){
        Tens <- samples[["U"]][s, ,]%*%samples[["R"]][s, v, ,]%*%t(samples[["V"]][s, ,])
        if(all(unique(val_df[["val_values"]], na.rm = T) %in% c(0, 1))){
          pROC::roc(val_df[["val_values"]], Tens[val_df[["val_indices"]]], auc = T) %>%
            magrittr::extract2("auc") %>%
            as.numeric()
        } else {
          sum( (Tens[val_df[["val_indices"]]] - val_df[["val_values"]])^2 ) %>%
            magrittr::divide_by(sum ((mean(val_df[["val_values"]]) - val_df[["val_values"]])^2)) %>%
            magrittr::subtract(1, .)
        }
      }) %>%
        data.frame(mean = mean(.),
                   sd = sd(.),
                   lwr = quantile(., 0.025),
                   upr = quantile(., 0.975)) %>%
        dplyr::select(mean, sd, lwr, upr) %>%
        dplyr::distinct()
    }) %>%
      magrittr::set_names(names(validation_set))
  } else if(type == "colwise"){
    plyr::alply(1:length(validation_set), 1, function(v){
      col_names <- colnames(tensor_list[[v]])
      val_df <- validation_set[[v]] %>%
                  dplyr::mutate(row_id = ((val_indices - 1)%%dim(samps[["U"]])[2]) + 1,
                                col_id = floor((val_indices - 1)/dim(samps[["U"]])[2]) + 1)
      plyr::aaply(1:n_samps, 1, function(s){
        Tens <- samples[["U"]][s, ,]%*%samples[["R"]][s, v, ,]%*%t(samples[["V"]][s, ,])
        plyr::aaply(unique(val_df[["col_id"]]), 1, function(cid){
          col_df <- val_df %>% dplyr::filter(col_id == cid)
          if(all(unique(col_df[["val_values"]], na.rm = T) %in% c(0, 1))){
            if(length(unique(col_df[["val_values"]], na.rm = T)) > 1){
              pROC::roc(col_df[["val_values"]], Tens[col_df[["val_indices"]]], auc = T) %>%
                magrittr::extract2("auc") %>%
                as.numeric()
            } else {
              NA
            }
          } else {
            sum( (Tens[col_df[["val_indices"]]] - col_df[["val_values"]])^2 ) %>%
              magrittr::divide_by(sum ((mean(tensor_list[[v]][setdiff(1:nrow(tensor_list[[v]]), col_df[["row_id"]]), cid], na.rm = T) - col_df[["val_values"]])^2)) %>%
              magrittr::subtract(1, .)
          }
        })
      }, .progress = "text") %>%
        plyr::adply(2, function(u){
          data.frame(mean = mean(u, na.rm = T),
                     sd = sd(u, na.rm = T),
                     lwr = quantile(u, 0.025, na.rm = T),
                     upr = quantile(u, 0.975, na.rm = T))
        }) %>%
        dplyr::mutate(col_name = col_names[unique(val_df[["col_id"]])]) %>%
        dplyr::select(col_name, mean, sd, lwr, upr) %>%
        dplyr::distinct()
    }) %>%
      magrittr::set_names(names(validation_set))
  }

}
