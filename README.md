# flexBCF
<!-- badges: start -->
  [![R-CMD-check](https://github.com/skdeshpande91/flexBCF/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/skdeshpande91/flexBCF/actions/workflows/R-CMD-check.yaml)
  [![CRAN status](https://www.r-pkg.org/badges/version/flexBCF)](https://CRAN.R-project.org/package=flexBCF)


<!-- badges: end -->


**flexBCF** is a faster and more flexible implementation of Bayesian Causal Forests.
This implementation was created as part of the 2022 American Causal Inference Competition (ACIC) Data Challenge.
Please see [our paper](https://doi.org/10.1353/obs.2023.0024) for more details.

Note that **flexBCF** implements a slightly different model than what is provided by default in the original **bcf** package (available at [this link](https://github.com/jaredsmurray/bcf)). Namely, it (i) does not use half-Normal or half-Cauchy priors for scale parameters; (ii) is not invariance to re-coding of the treatment indicator; and (iii) uses the same regression tree prior $\mu$ and $\tau$ (in contrast to **bcf** which more strongly regularizes $\tau$-trees towards shallow trees). If you require such functionality, please contact the package maintainer.

You can install **flexBCF** using
```
devtools::install_github(repo = "skdeshpande91/flexBCF")
```
