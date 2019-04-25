#'
#'
#'
subset_tensor <- function(tensor, idx1, idx2){
  
  require(plyr)
  require(magrittr)
  
  plyr::llply(tensor, function(m){
    
    m[idx1, idx2]
    
  }) %>%
    magrittr::set_names(names(tensor))
  
}
