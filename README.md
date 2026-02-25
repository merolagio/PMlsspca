
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Package PMlsspca

<!-- badges: start -->

[![R-CMD-check](https://github.com/merolagio/PMlsspca/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/merolagio/PMlsspca/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of PMlsspca is to help reproduce the Least Squares Sparse
Principal Components Analysis (LS SPCA) examples in the paper published
in …………… The main function is lsspca(). Methods are print, summery and
plot. Functions are new.spca (to create an \`spca’ object from a set of
loadings and aggregate_by_scale to visualize the contribution by scale.

## Installation

You can install the development version of PMlsspca from
[GitHub](https://github.com/) with:

``` r
# `remotes' is the lightest alternative
install.packages("remotes")
remotes::install_github("merolagio/PMlsspca")
#or
#install.packages("devtools")
devtools::install_github("merolagio/PMlsspca")
# or
# install.packages("pak")
pak::pak("merolagio/PMlsspca")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(PMlsspca)
#> Loading required package: Rcpp
#> Loading required package: RcppEigen
#> Loading required package: ggplot2
#> 
#> Attaching package: 'PMlsspca'
#> The following object is masked from 'package:stats':
#> 
#>     screeplot
## basic example code
#data(msscq)
#mylsspca = lsspca(ms, alpha = 0.95, ncomp = 4)
```

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->

<!-- ```{r cars} -->

<!-- summary(cars) -->

<!-- ``` -->

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. -->

<!-- You can also embed plots, for example: -->

<!-- ```{r pressure, echo = FALSE} -->

<!-- plot(pressure) -->

<!-- ``` -->

<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
