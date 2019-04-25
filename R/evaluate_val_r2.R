#'
#'
#'
evaluate_val_r2 <- function(validation_set, samples){
  #
  require(plyr)
  require(magrittr)
  #
  n_samps <- samples[["U"]] %>% nrow()
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
  
}