#'
#'
#'
make_validation_set <- function(tensor_list, val_pct = 0.05){
  #
  require(plyr)
  require(magrittr)
  #
  tensor_list %>%
    plyr::llply(function(m){

      numeric_indices <- which(!is.na(m))
      val_indices <- sample(numeric_indices, round(length(numeric_indices)*val_pct))
      data.frame(val_indices = val_indices,
                 val_values = m[val_indices])

    })

}

#'
#'
#'
hide_validation_set <- function(tensor_list, validation_set){
  #
  require(plyr)
  require(magrittr)
  #
  names(tensor_list) %>%
    plyr::alply(1, function(k){

      m <- tensor_list[[k]]
      m[validation_set[[k]][["val_indices"]]] <- NA
      m

    }) %>%
    magrittr::set_names(names(tensor_list))
}
