## Bayesian Tensor Factorization

Code for fitting GLRAM model with Gibbs sampling.

## Installation

To install, run the following command from the R console:

```{r}
devtools::install_github("j-g-b/BTF")
```

## Demo

First, define a function to simulate a three-way array of data with dimensions N, M, and K.

```{r}
simulate_tensor_data <- function(N, M, D, K){
  #
  U <- rnorm(N*D) %>%
        matrix(nrow = N)
  V <- rnorm(M*D) %>%
        matrix(nrow = M)
  R <- plyr::alply(1:K, 1, function(d){
        rnorm(D*D) %>% matrix(nrow = D)
       })
  SigmaSq <- 1/rgamma(K, 5/2, 1/2)
  #
  X <- plyr::llply(1:K, function(k){
        (U%*%R[[k]]%*%t(V)) %>%
          magrittr::set_rownames(paste0("CL", 1:N)) %>%
          magrittr::set_colnames(paste0("GN", 1:M)) %>%
          magrittr::add(matrix(rnorm(N*M, sd = sqrt(SigmaSq[k])), nrow = N))
       })
  #
  return(list(U = U, V = V, R = R, SigmaSq = SigmaSq, X = X))
}
```

Next choose dimensions and extract the components of the simulated three-way array X.

```{r}
#
require(plyr)
require(magrittr)
require(tidyverse)
#
N <- 100
M <- 500
D <- 5
K <- 5
#
dataset <- simulate_tensor_data(N, M, D, K)
#
U <- dataset[["U"]]
V <- dataset[["V"]]
R <- dataset[["R"]]
SigmaSq <- dataset[["SigmaSq"]]
X <- dataset[["X"]]
```

Next, run the tensor factorization to get samples of the embeddings U and V:

```{r}
#
Result <- BTF::run_btf(n_samples = 500, tensor_list = X, 
                       d1 = D, d2 = D, save_dir = "alt_samps", 
                       burn = 25, thin = 5)
#
Samples <- BTF::collect_samples("alt_samps/")
```