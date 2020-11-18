
# pedmod: Pedigree Models

TODO: write about the package.

## Example

First, we source a file to get a function to simulate a data set with a
maternal effect and a genetic effect:

``` r
# souce the file to get the simulation function
source(system.file("gen-pedigree-data.R", package = "pedmod"))

# simulate a data set
set.seed(68167102)
dat <- sim_pedigree_data(n_families = 1000L)
#> Loading required package: Matrix
#> Loading required package: quadprog

# distribution of family sizes
table(sapply(dat$sim_data, function(x) length(x$y)))
#> 
#>   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23 
#>  13  12   7   7   5   8  17  24  46 121 148 198 166  98  73  35  14   8

# total number of observations
sum(sapply(dat$sim_data, function(x) length(x$y)))
#> [1] 16665

# average event rate
mean(unlist(sapply(dat$sim_data, `[[`, "y")))
#> [1] 0.2354635

# TODO: show more about the families
```

Then we perform the model estimation:

``` r
# the true parameters are
dat$beta
#> (Intercept)          X1          X2 
#>        -1.0         0.3         0.2
dat$sc
#>   Gentic Maternal 
#>     0.50     0.33

# prepare the data to pass to the functions in the package
dat_arg <- lapply(dat$sim_data, function(x){
  # we need the following for each family: 
  #   y: the zero-one outcomes.
  #   X: the design matrix for the fixed effects. 
  #   scale_mats: list with the scale matrices for each type of effect.
  list(y = as.numeric(x$y), X = x$X,
       scale_mats = list(x$rel_mat, x$met_mat))
})

# create a C++ object
library(pedmod)
ll_terms <- get_pedigree_ll_terms(dat_arg, max_threads = 4L)

# get the starting values. This is very fast
y <- unlist(lapply(dat_arg, `[[`, "y"))
X <- do.call(rbind, lapply(dat_arg, `[[`, "X"))
start_fit <-  glm.fit(X, y, family = binomial("probit"))

# log-likelihood at the starting values without random effects
-sum(start_fit$deviance) / 2     
#> [1] -8825.879
(beta <- start_fit$coefficients) # starting values for fixed effects 
#> (Intercept)          X1             
#>  -0.7413368   0.2085147   0.1518448

# start at moderate sized scale parameters
sc <- rep(log(.2), 2)

# check log likelihood at the starting values. First we assign a function 
# to approximate the log likelihood and the gradient
fn <- function(par, seed = 1L, rel_eps = 1e-2, use_aprx = TRUE, 
               n_threads = 4L, indices = seq_along(dat_arg), 
               maxvls = 10000L){
  set.seed(seed)
  -eval_pedigree_ll(
    ll_terms, par = par, maxvls = maxvls, abs_eps = 0, rel_eps = rel_eps, 
    minvls = 1000L, use_aprx = use_aprx, n_threads = n_threads, 
    indices = indices)
}
gr <- function(par, seed = 1L, rel_eps = 1e-2, use_aprx = TRUE, 
               n_threads = 4L, indices = seq_along(dat_arg), 
               maxvls = 10000L){
  set.seed(seed)
  out <- -eval_pedigree_grad(
    ll_terms, par = par, maxvls = maxvls, abs_eps = 0, rel_eps = rel_eps, 
    minvls = 1000L, use_aprx = use_aprx, n_threads = n_threads, 
    indices = indices)
  structure(c(out), value = -attr(out, "logLik"))
}

# check output at the starting values
system.time(ll <- -fn(c(beta, sc)))
#>    user  system elapsed 
#>   0.966   0.000   0.244
ll # the log likelihood at the starting values
#> [1] -8687.09
system.time(gr_val <- gr(c(beta, sc)))
#>    user  system elapsed 
#>   3.106   0.003   0.809
gr_val # the gradient at the starting values
#> [1]  635.879395 -175.241684  -83.048473   17.660212   -8.328734
#> attr(,"value")
#> [1] 8687.083

# variance of the approximation
sd(sapply(1:25, function(seed) fn(c(beta, sc), seed = seed)))
#> [1] 0.03049123

# verify the gradient (may not be exactly equal due to MC error)
numDeriv::grad(fn, c(beta, sc))
#> [1]  636.312515 -175.447580  -83.008941   17.578706   -8.509126

# optimize the log likelihood approximation
system.time(opt <- optim(c(beta, sc), fn, gr, method = "BFGS"))
#>    user  system elapsed 
#> 222.061   0.000  56.587

# w/ higher precision
system.time(opt_prec <- optim(opt$par, fn, gr, method = "BFGS", 
                              rel_eps = 1e-3, maxvls = 25000L))
#>    user  system elapsed 
#> 453.707   0.000 117.270
```

The output from the optimization is shown below:

``` r
-opt$value      # the maximum log likelihood
#> [1] -8620.118
opt$convergence # check convergence
#> [1] 0

# w/ higher precision
-opt_prec$value
#> [1] -8620.088
opt_prec$convergence
#> [1] 0

# compare the estimated fixed effects with the true values
rbind(truth                     = dat$beta, 
      estimated                 = head(opt$par, length(dat$beta)), 
      `estimated (higher prec)` = head(opt_prec$par, length(dat$beta)))
#>                         (Intercept)        X1        X2
#> truth                     -1.000000 0.3000000 0.2000000
#> estimated                 -1.021337 0.2874056 0.2086335
#> estimated (higher prec)   -1.037563 0.2922601 0.2124689

# compare estimated scale parameters with the true values
rbind(truth                     = dat$sc, 
      estimated                 = exp(tail(opt$par     , length(dat$sc))), 
      `estimated (higher prec)` = exp(tail(opt_prec$par, length(dat$sc))))
#>                            Gentic  Maternal
#> truth                   0.5000000 0.3300000
#> estimated               0.6286649 0.2856720
#> estimated (higher prec) 0.7071673 0.2672134
```

### The Multivariate Normal CDF Approximation

We compare the multivariate normal CDF approximation in this section
with the approximation from the mvtnorm package. The same algorithm is
used but the version in this package is re-written in C++ and differs
slightly. Moreover, we have implemented an approximation of the standard
normal CDF and its inverse which reduces the computation time as we will
show below.

``` r
#####
# settings for the simulation study
library(mvtnorm)
library(pedmod)
library(microbenchmark)
set.seed(78459126)
n <- 5L         # number of variables to interate out
rel_eps <- 1e-4 # the relative error to use

#####
# run the simulation study
sim_res <- replicate(expr = {
  # simulate covariance matrix and the upper bound
  S <- drop(rWishart(1L, 2 * n, diag(n) / 2 / n))
  u <- rnorm(n)
  
  # function to use pmvnorm
  use_mvtnorm <- function(rel_eps)
    pmvnorm(upper = u, sigma = S, algorithm = GenzBretz(
      abseps = 0, releps = rel_eps, maxpts = 1e7))
  
  # function to use this package
  use_mvndst <- function(use_aprx = FALSE)
    mvndst(lower = rep(-Inf, n), upper = u, mu = rep(0, n), 
           sigma = S, use_aprx = use_aprx, abs_eps = 0, rel_eps = rel_eps,
           maxvls = 1e7)

  # get a very precise estimate
  truth <- use_mvtnorm(rel_eps / 100)
  
  # computes the error with repeated approximations and compute the time it
  # takes
  n_rep <- 5L
  run_n_time <- function(expr){
    expr <- substitute(expr)
    ti <- get_nanotime()
    res <- replicate(n_rep, eval(expr))
    ti <- get_nanotime() - ti
    err <- (res - truth) / truth
    c(SE = sqrt(sum(err^2) / n_rep), time = ti / n_rep / 1e9)
  }
  
  mvtnorm_res        <- run_n_time(use_mvtnorm(rel_eps))
  mvndst_no_aprx_res <- run_n_time(use_mvndst(FALSE))
  mvndst_w_aprx_res  <- run_n_time(use_mvndst(TRUE))
  
  # return 
  rbind(mvtnorm            = mvtnorm_res, 
        `mvndst (no aprx)` = mvndst_no_aprx_res, 
        `mvndst (w/ aprx)` = mvndst_w_aprx_res)
}, n = 100, simplify = "array")
```

They have about the same average relative error as expected:

``` r
rowMeans(sim_res[, "SE", ])
#>          mvtnorm mvndst (no aprx) mvndst (w/ aprx) 
#>     0.0000282989     0.0000316296     0.0164280867
par(mar = c(5, 4, 1, 1))
boxplot(t(sim_res[, "SE", ]))
```

<img src="man/figures/README-show_averge_rel_err-1.png" width="100%" />

``` r

# without extreme values
sum(keep <- colSums(sim_res[, "SE", ] > .5) < 1)
#> [1] 99
rowMeans(sim_res [, "SE", keep])
#>          mvtnorm mvndst (no aprx) mvndst (w/ aprx) 
#>     2.831528e-05     3.176110e-05     3.157978e-05
boxplot(t(sim_res[, "SE", keep]))
```

<img src="man/figures/README-show_averge_rel_err-2.png" width="100%" />

The new implementation is faster when the approximation is used:

``` r
rowMeans(sim_res[, "time", ])
#>          mvtnorm mvndst (no aprx) mvndst (w/ aprx) 
#>       0.01707832       0.01653099       0.01094785
par(mar = c(5, 4, 1, 1))
boxplot(t(sim_res[, "time", ]))
```

<img src="man/figures/README-use_new_impl-1.png" width="100%" />
