#'
#'
#'
nearest_neighbors <- function(U, k){

  require(FNN)
  require(magrittr)
  require(plyr)

  FNN::get.knn(U, k = k) %>%
    magrittr::extract2("nn.index") %>%
    magrittr::set_rownames(row.names(U)) %>%
    plyr::aaply(1, function(x){
      magrittr::extract(row.names(U), x)
    })

}

query_distance <- function(query, U, row_names, summary = T){
  require(magrittr)

  if(length(dim(U)) < 3){
    U %<>% array(c(1, dim(U)))
  }

  plyr::alply(query, 1, function(q){
    all_res <- plyr::aaply(U, 1, function(u){

      u %<>% magrittr::set_rownames(row_names)
      u %>%
        apply(1, function(uu){
          sum(uu*u[q, ]) / (sqrt(sum(uu^2))*sqrt(sum(u[q, ]^2)))
        })
    })

    if(length(dim(all_res)) < 2){
      all_res %<>% array(c(1, length(all_res))) %>% magrittr::set_colnames(row_names)
    }

    if(summary){
      all_res %>%
        plyr::adply(2, function(col){
          data.frame(mean = mean(col),
                     sd = sd(col),
                     lwr = quantile(col, 0.025),
                     upr = quantile(col, 0.975))
        }) %>%
        dplyr::mutate(query = q) %>%
        dplyr::rename(reference = X1) %>%
        dplyr::select(query, reference, mean, sd, lwr, upr) %>%
        dplyr::arrange(mean)
    } else {
      all_res
    }
  }) %>%
    magrittr::set_names(query)

}

nearest_neighbor_distribution <- function(U, k, row_names){

  require(FNN)
  require(magrittr)
  require(plyr)

  interim_res <- plyr::aaply(1:nrow(U), 1, function(i){

    FNN::get.knn(U[i,,], k = k) %>%
      magrittr::extract2("nn.index") %>%
      magrittr::set_rownames(row_names) %>%
      plyr::aaply(1, function(x){
        magrittr::extract(row_names, x)
      })

  })

  plyr::adply(1:(dim(interim_res)[2]), 1, function(i){
      m <- interim_res[,i,]
      plyr::adply(1:ncol(m), 1, function(j){
        plyr::count(m[,j]) %>%
          dplyr::mutate(neighbor_k = j) %>%
          dplyr::mutate(anchor = row_names[i]) %>%
          dplyr::mutate(freq = freq / nrow(m))
      })
    }) %>%
    dplyr::select(anchor, x, neighbor_k, freq) %>%
    plyr::arrange(anchor, neighbor_k, 1/freq)

}
