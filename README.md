RZiMM-scRNA
================

**RZiMM-scRNA**, a regularized zero-inflated mixture model designed for
sparse single-cell RNA-seq data with sample heterogeneity (Paper to be
submitted).

RZiMM-scRNA take within-sample correlation into consideration,
simultaneously cluster the cells into groups and identify differentially
expressed genes among groups.

# Installation

In R,

``` r
devtools::install_github("SkadiEye/RZiMM")
```

# Example (Normalized)

``` r
suppressMessages(library(RZiMM))
## Load example data
data("sim_x")
x <- sim_x$x
group_ind <- sim_x$cluster
sample_index <- sim_x$sample_index
imp_v <- sim_x$imp_gene_list

set.seed(4321)
## Run RZiMM-scRNA (adjusting for sample correlation)
mod_rzimm_sc <- RZiMM(x, sample_index, n_group = 4, lambda1 = 30)
## Run RZiMM-Naive (not adjusting for sample correlation)
mod_rzimm_na <- RZiMM(x, n_group = 4, lambda1 = 30)
## Run RZiMM-GMM (Gaussian-Mixture)
mod_rzimm_gm <- RZiMM_GMM(x, sample_index, n_group = 4, lambda1 = 0)

## Clustering performance
rbind(
  data.frame(
    method = "RZiMM-scRNA", 
    clustering_scores(group_ind, mod_rzimm_sc@cluster), 
    AUC = perf_stat(mod_rzimm_sc@importance, 0.1, imp_v)$AUC), 
  data.frame(
    method = "RZiMM-Naive", 
    clustering_scores(group_ind, mod_rzimm_na@cluster), 
    AUC = perf_stat(mod_rzimm_na@importance, 0.1, imp_v)$AUC), 
  data.frame(
    method = "RZiMM-GMM", 
    clustering_scores(group_ind, mod_rzimm_gm@cluster), 
    AUC = perf_stat(mod_rzimm_gm@importance, 1, imp_v)$AUC))
```

    ##        method       ARI       NMI Homogeneity       AUC
    ## 1 RZiMM-scRNA 1.0000000 1.0000000   1.0000000 1.0000000
    ## 2 RZiMM-Naive 0.5945923 0.7081437   0.7380899 0.9994222
    ## 3   RZiMM-GMM 1.0000000 1.0000000   1.0000000 1.0000000

# Example (Count)

``` r
data("sim_y")
x <- sim_y$x
group_ind <- sim_y$cluster
sample_index <- sim_y$sample_index
imp_v <- sim_y$imp_gene_list

set.seed(10086)
## Run RZiMM-NB 
mod_rzimm_nb <- RZiMM_Count(x, sample_index, n_group = 4, dist = "zinb")
## Run RZiMM-Pois 
mod_rzimm_ps <- RZiMM_Count(x, sample_index, n_group = 4, dist = "zip")
## Run RZiMM-scRNA 
mod_rzimm_sc <- RZiMM(x, sample_index, n_group = 4)

## Clustering performance
rbind(
  data.frame(
    method = "RZiMM-NB", 
    clustering_scores(group_ind, mod_rzimm_nb@cluster), 
    AUC = perf_stat(mod_rzimm_nb@importance, 0.1, imp_v)$AUC), 
  data.frame(
    method = "RZiMM-Pois", 
    clustering_scores(group_ind, mod_rzimm_ps@cluster), 
    AUC = perf_stat(mod_rzimm_ps@importance, 0.1, imp_v)$AUC), 
  data.frame(
    method = "RZiMM-scRNA", 
    clustering_scores(group_ind, mod_rzimm_sc@cluster), 
    AUC = perf_stat(mod_rzimm_sc@importance, 1, imp_v)$AUC))
```

    ##        method       ARI       NMI Homogeneity AUC
    ## 1    RZiMM-NB 0.9663742 0.9419806   0.9419342   1
    ## 2  RZiMM-Pois 0.9711466 0.9522562   0.9521428   1
    ## 3 RZiMM-scRNA 0.9661785 0.9471121   0.9467537   1

# References
