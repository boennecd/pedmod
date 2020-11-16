## ---- include = FALSE---------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)


## -----------------------------------------------------------------------------------
n <- 5L
set.seed(1)
S <- drop(rWishart(1L, 2 * n, diag(n) / 2 / n))
u <- rnorm(n)

library(mvtnorm)
library(pedmod)
pmvnorm(upper = u, sigma = S)
mvndst(lower = rep(-Inf, n), upper = u, mu = rep(0, n), 
       sigma = S)

