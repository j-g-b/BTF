## Bayesian Tensor Factorization

Code for Bayesian inference in the tensor probability model described by ["Smaller $p$-values in genomics studies using distilled historical information"](https://arxiv.org/abs/2004.07887).

## Installation

To install, run the following command from the R console:

```{r}
devtools::install_github("j-g-b/BTF")
```

## Running the simulation studies

The scripts in [scripts/simulation_studies](https://github.com/j-g-b/BTF/tree/master/scripts/simulation_studies) should run without the need for additional data files. To save the resulting figures, make sure to create a directory called `figs` before running.

## Running the genomics data examples

### (1) Download the samples

From this [figshare link](https://figshare.com/projects/Smaller_p-values_in_genomics_studies_using_distilled_historical_information/79287) download the .zip file in the dataset called "Posterior samples." Unzip the file to extract all of the .csv files and put these in a separate directory. For example these could be stored in a directory called `"samps"`.

### (2) Choose a script

Select an example dataset from the scripts in this repo under [scripts/genomics_examples](https://github.com/j-g-b/BTF/tree/master/scripts/genomics_examples). Then follow the [figshare link](https://figshare.com/projects/Smaller_p-values_in_genomics_studies_using_distilled_historical_information/79287) and download the corresponding data file. For example, if running `fabp_corsello_et_al.R`, download the "Corsello et al. example" dataset. Place the script in the same directory as the downloaded dataset.

### (3) Modify the `samps_dir`

At the top of the example script, there is a line that sets the name of the directory where the posterior samples are saved. Modify the example script to point to the place where the posterior samples were saved in step (1).

For example, the top lines of the `fabp_kory_et_al.R` script look like this

```{r}
#
samps_dir <- "../../alt_samps"
#
library(tidyverse)
library(magrittr)
#
K562_table <- readr::read_csv("K562 SHMT1-null.csv") %>%
                dplyr::filter(!grepl("INTERGENIC", sgRNA)) %>%
                dplyr::mutate(gene = gsub("_.*", "", gsub("sg", "", sgRNA)),
                              replicate = as.integer(gsub(".*_", "", sgRNA)),
                              full_lfc = log2(((`full media` + 1)/sum(`full media` + 1, na.rm = T)) / ((`initial` + 1)/sum(`initial` + 1, na.rm = T))),
                              minus_lfc = log2(((`minus serine` + 1)/sum(`minus serine` + 1, na.rm = T)) / ((`initial` + 1)/sum(`initial` + 1, na.rm = T)))) %>%
                dplyr::select(-sgRNA) %>%
                dplyr::group_by(gene) %>%
                dplyr::summarise(r = sum(!is.na(minus_lfc - full_lfc)),
                                 s = sd(minus_lfc - full_lfc, na.rm = T),
                                 m = mean(minus_lfc - full_lfc, na.rm = T)) %>%
                dplyr::filter(r > 1)
```

Simply change the name of the directory to stored in `samps_dir` by modifying that line to read, for example

```{r}
samps_dir <- "my_samps"
```

### (4) Run the script

The example script is ready to run. Make sure to create a directory called `figs` in order to save the resulting figure properly.

## Running Gibbs sampling on the Cancer Dependency Map data

### (1) Download the data from DepMap

Visit the [depmap.org downloads page](https://depmap.org/portal/download/) to find the gene expression, mutation, CRISPR, and RNAi datasets. Download the following files:

1. `CCLE_expression.csv`
2. `CCLE_mutations.csv`
3. `Achilles_gene_effect.csv`
4. `D2_combined_gene_dep_scores.csv`
5. `sample_info.csv`

and put them into a directory called `data`.

### (2) Download the gene groups from figshare

Download the files in the "Gene groups" dataset found at this [figshare project](https://figshare.com/projects/Smaller_p-values_in_genomics_studies_using_distilled_historical_information/79287) and place them in the `data` directory.

### (3) Prepare the data

Run the script [`PrepareInputs.R`](https://github.com/j-g-b/BTF/blob/master/scripts/PrepareInputs.R). This will create `.rds` objects of the DepMap data files.

### (4) Run Gibbs sampling

Run the script [`generate_posterior_samples.R`](https://github.com/j-g-b/BTF/blob/master/scripts/generate_posterior_samples.R), perhaps first modifying the last command to change the name of the directory in which posterior samples are stored:

```{r}
X <- BTF::run_btf(n_samples = S, tensor_list = TensorList,
                  d1 = d1, d2 = d2, save_dir = "my_choice_of_directory_name",
                  burn = 1000, thin = 40, matrix_type = c(2, 0, 0, 1))
```

## Demo for running the Gibbs sampling algorithm on generic three-way tensor data

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
                       d1 = D, d2 = D, save_dir = "samps", 
                       burn = 25, thin = 5)
#
Samples <- BTF::collect_samples("samps/")
```

## Getting posterior means from posterior samples

Due to rotation and scale ambiguity, the posterior samples of the model parameters need to be aligned before taking posterior summaries.

Suppose that the name of the directory where the posterior sample files have been written is stored in the variable `samps_dir`. Then matrices `U` and `V` containing the posterior means of the row and column features can be obtained using the following commands:

```{r}
samps <- BTF::collect_samples(samps_dir, last = 100) %>% BTF::align_samples()
V <- apply(samps$V, c(2,3), mean)
U <- apply(samps$U, c(2,3), mean)
```
