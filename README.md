
# pedmod: Pedigree Models

[![R-CMD-check](https://github.com/boennecd/pedmod/workflows/R-CMD-check/badge.svg)](https://github.com/boennecd/pedmod/actions)

The pedmod package provides functions to estimate models for pedigree
data. Particularly, the package provides functions to estimate mixed
models of the form:

  
![\\begin{align\*}&#10;Y\_{ij} \\mid \\epsilon\_{ij} = e &#10; &\\sim
\\text{Bin}(\\Phi(\\vec\\beta^\\top\\vec x\_{ij} + e), 1)
\\\\&#10;\\vec\\epsilon\_i = (\\epsilon\_{i1}, \\dots,
\\epsilon\_{in\_i})^\\top &\\sim&#10; N^{(n\_i)}\\left(\\vec 0,
\\sum\_{l = 1}^K\\sigma\_l^2 C\_{il}&#10;
\\right)&#10;\\end{align\*}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign%2A%7D%0AY_%7Bij%7D%20%5Cmid%20%5Cepsilon_%7Bij%7D%20%3D%20e%20%0A%20%20%26%5Csim%20%5Ctext%7BBin%7D%28%5CPhi%28%5Cvec%5Cbeta%5E%5Ctop%5Cvec%20x_%7Bij%7D%20%2B%20e%29%2C%201%29%20%5C%5C%0A%5Cvec%5Cepsilon_i%20%3D%20%28%5Cepsilon_%7Bi1%7D%2C%20%5Cdots%2C%20%5Cepsilon_%7Bin_i%7D%29%5E%5Ctop%20%26%5Csim%0A%20%20N%5E%7B%28n_i%29%7D%5Cleft%28%5Cvec%200%2C%20%5Csum_%7Bl%20%3D%201%7D%5EK%5Csigma_l%5E2%20C_%7Bil%7D%0A%20%20%5Cright%29%0A%5Cend%7Balign%2A%7D
"\\begin{align*}
Y_{ij} \\mid \\epsilon_{ij} = e 
  &\\sim \\text{Bin}(\\Phi(\\vec\\beta^\\top\\vec x_{ij} + e), 1) \\\\
\\vec\\epsilon_i = (\\epsilon_{i1}, \\dots, \\epsilon_{in_i})^\\top &\\sim
  N^{(n_i)}\\left(\\vec 0, \\sum_{l = 1}^K\\sigma_l^2 C_{il}
  \\right)
\\end{align*}")  

where
![Y\_{ij}](https://render.githubusercontent.com/render/math?math=Y_%7Bij%7D
"Y_{ij}") is the binary outcome of interest for individual
![j](https://render.githubusercontent.com/render/math?math=j "j") in
family/cluster
![i](https://render.githubusercontent.com/render/math?math=i "i"),
![\\vec
x\_{ij}](https://render.githubusercontent.com/render/math?math=%5Cvec%20x_%7Bij%7D
"\\vec x_{ij}") is the individual’s known covariates,
![\\Phi](https://render.githubusercontent.com/render/math?math=%5CPhi
"\\Phi") is the standard normal distribution’s CDF, and
![\\text{Bin}](https://render.githubusercontent.com/render/math?math=%5Ctext%7BBin%7D
"\\text{Bin}") implies a binomial distribution such if ![z\\sim
\\text{Bin}(p,
n)](https://render.githubusercontent.com/render/math?math=z%5Csim%20%5Ctext%7BBin%7D%28p%2C%20n%29
"z\\sim \\text{Bin}(p, n)") then the density of
![z](https://render.githubusercontent.com/render/math?math=z "z") is:

  
![f(z) = \\begin{pmatrix} n \\\\ z
\\end{pmatrix}p^zp^{n-z}](https://render.githubusercontent.com/render/math?math=f%28z%29%20%3D%20%5Cbegin%7Bpmatrix%7D%20n%20%5C%5C%20z%20%5Cend%7Bpmatrix%7Dp%5Ezp%5E%7Bn-z%7D
"f(z) = \\begin{pmatrix} n \\\\ z \\end{pmatrix}p^zp^{n-z}")  

A different and equivalent way of writing the model is as:

  
![\\begin{align\*}&#10;Y\_{ij} \\mid \\epsilon\_{ij} = e &#10; &=
\\begin{cases}&#10; 1 & \\vec\\beta^\\top\\vec x\_{ij} + e \> 0
\\\\&#10; 0 & \\text{otherwise}&#10; \\end{cases}
\\\\&#10;\\vec\\epsilon\_i = (\\epsilon\_{i1}, \\dots,
\\epsilon\_{in\_i})^\\top &\\sim&#10; N^{(n\_i)}\\left(\\vec 0,
I\_{n\_i} + \\sum\_{l = 1}^K\\sigma\_l^2 C\_{il}&#10;
\\right)&#10;\\end{align\*}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign%2A%7D%0AY_%7Bij%7D%20%5Cmid%20%5Cepsilon_%7Bij%7D%20%3D%20e%20%0A%20%20%26%3D%20%5Cbegin%7Bcases%7D%0A%20%20%20%201%20%26%20%5Cvec%5Cbeta%5E%5Ctop%5Cvec%20x_%7Bij%7D%20%2B%20e%20%3E%200%20%5C%5C%0A%20%20%20%200%20%26%20%5Ctext%7Botherwise%7D%0A%20%20%20%20%5Cend%7Bcases%7D%20%5C%5C%0A%5Cvec%5Cepsilon_i%20%3D%20%28%5Cepsilon_%7Bi1%7D%2C%20%5Cdots%2C%20%5Cepsilon_%7Bin_i%7D%29%5E%5Ctop%20%26%5Csim%0A%20%20N%5E%7B%28n_i%29%7D%5Cleft%28%5Cvec%200%2C%20I_%7Bn_i%7D%20%2B%20%5Csum_%7Bl%20%3D%201%7D%5EK%5Csigma_l%5E2%20C_%7Bil%7D%0A%20%20%5Cright%29%0A%5Cend%7Balign%2A%7D
"\\begin{align*}
Y_{ij} \\mid \\epsilon_{ij} = e 
  &= \\begin{cases}
    1 & \\vec\\beta^\\top\\vec x_{ij} + e \> 0 \\\\
    0 & \\text{otherwise}
    \\end{cases} \\\\
\\vec\\epsilon_i = (\\epsilon_{i1}, \\dots, \\epsilon_{in_i})^\\top &\\sim
  N^{(n_i)}\\left(\\vec 0, I_{n_i} + \\sum_{l = 1}^K\\sigma_l^2 C_{il}
  \\right)
\\end{align*}")  

where
![I\_{n\_i}](https://render.githubusercontent.com/render/math?math=I_%7Bn_i%7D
"I_{n_i}") is the
![n\_i](https://render.githubusercontent.com/render/math?math=n_i "n_i")
dimensional identity matrix which comes from the unshared/individual
specific random effect. This effect is always included.

The
![C\_{il}](https://render.githubusercontent.com/render/math?math=C_%7Bil%7D
"C_{il}")s are known scale/correlation matrices where each of the
![l](https://render.githubusercontent.com/render/math?math=l "l")’th
types correspond to a type of effect. An arbitrary number of such
matrices can be passed to include e.g. a genetic effect, a maternal
effect, a paternal, an effect of a shared adult environment etc. A
typical example is that
![C\_{il}](https://render.githubusercontent.com/render/math?math=C_%7Bil%7D
"C_{il}") is two times the kinship matrix in which case we call:

  
![\\frac{\\sigma\_l^2}{1 + \\sum\_{k
= 1}^K\\sigma\_k^2}](https://render.githubusercontent.com/render/math?math=%5Cfrac%7B%5Csigma_l%5E2%7D%7B1%20%2B%20%5Csum_%7Bk%20%3D%201%7D%5EK%5Csigma_k%5E2%7D
"\\frac{\\sigma_l^2}{1 + \\sum_{k = 1}^K\\sigma_k^2}")  

the heritability (the proportion of the variance attributable to the
direct genetic effect). The scale parameters, the
![\\sigma\_k^2](https://render.githubusercontent.com/render/math?math=%5Csigma_k%5E2
"\\sigma_k^2")s, may be the primary interest in an analysis. We fix the
unshared random effect’ variance to be one which gives the one in the
denominator. Other authors prefer that the denominator is restricted to
be one. Both restriction can be applied to ensure identification and
they do not have any implications as the scale is arbitrary.

This package provides randomized quasi-Monte Carlo methods to
approximate the log marginal likelihood for these types of models with
an arbitrary number scale matrices,
![K](https://render.githubusercontent.com/render/math?math=K "K"), and
the derivatives with respect to
![(\\vec\\beta^\\top, 2\\log\\sigma\_1,\\dots, 2\\log\\sigma\_K)^\\top](https://render.githubusercontent.com/render/math?math=%28%5Cvec%5Cbeta%5E%5Ctop%2C%202%5Clog%5Csigma_1%2C%5Cdots%2C%202%5Clog%5Csigma_K%29%5E%5Ctop
"(\\vec\\beta^\\top, 2\\log\\sigma_1,\\dots, 2\\log\\sigma_K)^\\top")
(that is, we work with ![\\theta\_k
= 2\\log\\sigma\_k](https://render.githubusercontent.com/render/math?math=%5Ctheta_k%20%3D%202%5Clog%5Csigma_k
"\\theta_k = 2\\log\\sigma_k")). We have re-written the Fortran code by
Genz and Bretz (2002) in C++, made it easy to extend from a log marginal
likelihood approximation to other approximations such as the
derivatives, and added less precise but faster approximations of the
![\\Phi](https://render.githubusercontent.com/render/math?math=%5CPhi
"\\Phi") and
![\\Phi^{-1}](https://render.githubusercontent.com/render/math?math=%5CPhi%5E%7B-1%7D
"\\Phi^{-1}"). Our own experience suggests that using the latter has a
small effect on the precision of the result but can yield substantial
reduction in computation times for moderate sized families/clusters.

The approximation by Genz and Bretz (2002) have already been used to
estimate these types of models (Pawitan et al. 2004). However, not
having the gradients may slow down estimation substantially. Moreover,
our implementation supports computation in parallel which is a major
advantage given the availability of multi-core processors.

Since the implementation is easy to extend, possible extensions are:

1.  Survival times using mixed generalized survival models (Liu,
    Pawitan, and Clements 2017) with a similar random effect structure
    as the model shown above. This way, one avoids dichotomizing
    outcomes and can account for censoring.
2.  Generalized linear mixed model with binary, binomial, ordinal, or
    multinomial outcomes with a probit link. The method we use here may
    be beneficial if the number of random effects per cluster is not
    much smaller then the number observations in each cluster.

## Installation

The package can be installed from Github by calling:

``` r
remotes::install_github("boennecd/pedmod")
```

## Example

We start with a simple example only with a direct genetic effect. We
have one type of family which consists of two couples which are related
through one of the parents being cousins. The family is shown below.

``` r
# create the family we will use
fam <- data.frame(id = 1:10, sex = rep(1:2, 5L),
                  father = c(NA, NA, 1L, NA, 1L, NA, 3L, 3L, 5L, 5L), 
                  mother = c(NA, NA, 2L, NA, 2L, NA, 4L, 4L, 6L, 6L))

# plot the pedigree
library(kinship2)
ped <- with(fam, pedigree(id = id, dadid = father, momid = mother, sex = sex))
plot(ped)
```

<img src="man/figures/README-setup_simple-1.png" width="100%" />

We set the scale matrix to be two times the kinship matrix to model the
direct genetic effect. Each individual also have a standard normally
distributed covariate and a binary covariate. Thus, we can simulate a
data set with a function like:

``` r
# simulates a data set. 
# 
# Args:
#   n_fams: number of families.
#   beta: the fixed effect coefficients.
#   sig_sq: the scale parameter.
sim_dat <- function(n_fams, beta = c(-3, 1, 2), sig_sq = 3){
  # setup before the simulations
  Cmat <- 2 * kinship(ped)
  n_obs <- NROW(fam)
  Sig <- diag(n_obs) + sig_sq * Cmat
  Sig_chol <- chol(Sig)
  
  # simulate the data
  out <- replicate(
    n_fams, {
      # simulate covariates
      X <- cbind(`(Intercept)` = 1, Continuous = rnorm(n_obs), 
                 Binary = runif(n_obs) > .5)
      
      # assign the linear predictor + noise
      eta <- drop(X %*% beta) + drop(rnorm(n_obs) %*% Sig_chol)
      
      # return the list in the format needed for the package
      list(y = as.numeric(eta > 0), X = X, scale_mats = list(Cmat))
    }, simplify = FALSE)
  
  # add attributes with the true values and return 
  attributes(out) <- list(beta = beta, sig_sq = sig_sq)
  out
}
```

The model is

  
![\\begin{align\*}&#10; Y\_{ij} &= \\begin{cases} 1 & \\beta\_0 +
\\beta\_1 X\_{ij} + \\beta\_2 B\_{ij} + G\_{ij} + R\_{ij} \> 0 \\\\ 0 &
\\text{otherwise} \\end{cases} \\\\&#10; X\_{ij} &\\sim N(0, 1)
\\\\&#10; B\_{ij} &\\sim \\text{Bin}(0.5, 1) \\\\&#10; (G\_{i1}, \\dots
G\_{in\_{i}})^\\top &\\sim N^{(n\_i)}(\\vec 0, \\sigma^2 C\_{i1})
\\\\&#10; R\_{ij} &\\sim
N(0, 1)\\end{align\*}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign%2A%7D%0A%20Y_%7Bij%7D%20%26%3D%20%5Cbegin%7Bcases%7D%201%20%26%20%5Cbeta_0%20%2B%20%5Cbeta_1%20X_%7Bij%7D%20%2B%20%5Cbeta_2%20B_%7Bij%7D%20%2B%20G_%7Bij%7D%20%2B%20R_%7Bij%7D%20%3E%200%20%5C%5C%200%20%26%20%5Ctext%7Botherwise%7D%20%5Cend%7Bcases%7D%20%5C%5C%0A%20X_%7Bij%7D%20%26%5Csim%20N%280%2C%201%29%20%5C%5C%0A%20B_%7Bij%7D%20%26%5Csim%20%5Ctext%7BBin%7D%280.5%2C%201%29%20%5C%5C%0A%20%28G_%7Bi1%7D%2C%20%5Cdots%20G_%7Bin_%7Bi%7D%7D%29%5E%5Ctop%20%26%5Csim%20N%5E%7B%28n_i%29%7D%28%5Cvec%200%2C%20%5Csigma%5E2%20C_%7Bi1%7D%29%20%5C%5C%0A%20R_%7Bij%7D%20%26%5Csim%20N%280%2C%201%29%5Cend%7Balign%2A%7D
"\\begin{align*}
 Y_{ij} &= \\begin{cases} 1 & \\beta_0 + \\beta_1 X_{ij} + \\beta_2 B_{ij} + G_{ij} + R_{ij} \> 0 \\\\ 0 & \\text{otherwise} \\end{cases} \\\\
 X_{ij} &\\sim N(0, 1) \\\\
 B_{ij} &\\sim \\text{Bin}(0.5, 1) \\\\
 (G_{i1}, \\dots G_{in_{i}})^\\top &\\sim N^{(n_i)}(\\vec 0, \\sigma^2 C_{i1}) \\\\
 R_{ij} &\\sim N(0, 1)\\end{align*}")  

where
![C\_{i1}](https://render.githubusercontent.com/render/math?math=C_%7Bi1%7D
"C_{i1}") is two times the kinship matrix and
![X\_{ij}](https://render.githubusercontent.com/render/math?math=X_%7Bij%7D
"X_{ij}") and
![B\_{ij}](https://render.githubusercontent.com/render/math?math=B_%7Bij%7D
"B_{ij}") are observed covariates. We can now estimate the model with a
simulated data set as follows:

``` r
# simulate a data set
set.seed(27107390)
dat <- sim_dat(n_fams = 400L)

# perform the optimization. We start with finding the starting values
library(pedmod)
ll_terms <- get_pedigree_ll_terms(dat, max_threads = 4L)
system.time(start <- pedmod_start(ptr = ll_terms, data = dat, n_threads = 4L))
#>    user  system elapsed 
#>   1.774   0.017   0.469

# log-likelihood without the random effects and at the starting values
start$logLik_no_rng
#> [1] -1690
start$logLik_est # this is unreliably/imprecise
#> [1] -1619

# estimate the model
system.time(
  opt_out <- pedmod_opt(
    ptr = ll_terms, par = start$par, abs_eps = 0, use_aprx = TRUE, 
    n_threads = 4L, 
    maxvls = 25000L, rel_eps = 1e-3, minvls = 5000L))
#>    user  system elapsed 
#> 102.270   0.024  25.647
```

The results are shown below:

``` r
# parameter estimates versus the truth
rbind(opt_out       = head(opt_out$par, -1), 
      opt_out_quick = head(start  $par, -1), 
      truth         = attr(dat, "beta"))
#>               (Intercept) Continuous Binary
#> opt_out            -2.880     0.9714  1.883
#> opt_out_quick      -2.868     0.9944  1.873
#> truth              -3.000     1.0000  2.000
c(opt_out       = exp(tail(opt_out$par, 1)), 
  opt_out_quick = exp(tail(start  $par, 1)), 
  truth         = attr(dat, "sig_sq"))
#>       opt_out opt_out_quick         truth 
#>         2.929         2.877         3.000

# log marginal likelihoods
print(start   $logLik_est, digits = 8) # this is unreliably/imprecise
#> [1] -1618.5072
print(-opt_out$value     , digits = 8)
#> [1] -1618.4042
```

We can compute a profile likelihood curve like this:

``` r
# assign the scale parameter at which we will evaluate the profile likelihood
rg <- range(exp(tail(opt_out$par, 1) / 2) * c(.5, 2),
            sqrt(attr(dat, "sig_sq")) * c(.9, 1.1))
sigs <- seq(rg[1], rg[2], length.out = 10)
sigs <- sort(c(sigs, exp(tail(opt_out$par, 1) / 2)))

# compute the profile likelihood
ll_terms <- get_pedigree_ll_terms(dat, max_threads = 4L)
pl_curve_res <- lapply(sigs, function(sig){
  # set the parameters to pass
  beta <- start$beta_no_rng
  sig_sq_log <- 2 * log(sig)
  beta_scaled <- beta * sqrt(1 + sig^2)
  
  # optimize like before but using the fix argument
  opt_out_quick <- pedmod_opt(
    ptr = ll_terms, par = c(beta_scaled, sig_sq_log), maxvls = 1000L, 
    abs_eps = 0, rel_eps = 1e-2, minvls = 100L, use_aprx = TRUE, n_threads = 4L, 
    fix = length(beta) + 1L)
  opt_out <- pedmod_opt(
    ptr = ll_terms, par = c(opt_out_quick$par, sig_sq_log), abs_eps = 0, 
    use_aprx = TRUE, n_threads = 4L, fix = length(beta) + 1L,
    # we changed the parameters
    maxvls = 25000L, rel_eps = 1e-3, minvls = 5000L)
  
  # report to console and return
  message(sprintf("\nLog-likelihood %.5f (%.5f). Estimated parameters:", 
                  -opt_out$value, -opt_out_quick$value))
  message(paste0(capture.output(print(
    c(opt_out$par, Scale = sig))), collapse = "\n"))
  
  list(opt_out_quick = opt_out_quick, opt_out = opt_out)
})
```

We can construct an approximate 95% confidence interval using an
estimated cubic smoothing spline for the profile likelihood (more `sigs`
points may be needed to get a good estimate of the cubic smoothing
spline):

``` r
# check the critical values
alpha <- .05
crit_val <- qchisq(1 - alpha, 1)

# fit the cubic smoothing spline
pls <- -sapply(pl_curve_res, function(x) x$opt_out$value)
smooth_est <- smooth.spline(sigs, pls)

# check that we have values within the bounds
max_ml <- -opt_out$value
ll_diffs <- 2 * (max_ml - pls)
stopifnot(any(head(ll_diffs, length(ll_diffs) / 2) > crit_val), 
          any(tail(ll_diffs, length(ll_diffs) / 2) > crit_val))

# find the values
max_par <- tail(opt_out$par, 1)
lb <- uniroot(function(x) 2 * (max_ml - predict(smooth_est, x)$y) - crit_val, 
              c(min(sigs)       , exp(max_par / 2)))$root
ub <- uniroot(function(x) 2 * (max_ml - predict(smooth_est, x)$y) - crit_val, 
              c(exp(max_par / 2), max(sigs)))$root

# the confidence interval 
c(lb, ub)
#> [1] 1.259 2.528
c(lb, ub)^2 # on the variance scale
#> [1] 1.586 6.391
```

A caveat is that issues with the
![\\chi^2](https://render.githubusercontent.com/render/math?math=%5Cchi%5E2
"\\chi^2") approximation may arise on the boundary of the scale
parameter (![\\sigma
= 0](https://render.githubusercontent.com/render/math?math=%5Csigma%20%3D%200
"\\sigma = 0"); see <https://stats.stackexchange.com/a/4894/81865>).
Notice that the above may fail if the estimated profile likelihood is
not smooth e.g. because of convergence issues. We can plot the profile
likelihood and highlight the critical value as follows:

``` r
par(mar = c(5, 5, 1, 1))
plot(sigs, pls, bty = "l",
     pch = 16, xlab = expression(sigma), ylab = "Profile likelihood")
grid()
lines(predict(smooth_est, seq(min(sigs), max(sigs), length.out = 100)))
abline(v = exp(tail(opt_out$par, 1) / 2), lty = 2) # the estimate
abline(v = sqrt(attr(dat, "sig_sq")), lty = 3) # the true value
abline(v = lb, lty = 3) # mark the lower bound
abline(v = ub, lty = 3) # mark the upper bound
abline(h = max_ml - crit_val / 2, lty = 3) # mark the critical value
```

<img src="man/figures/README-plot_simple_ex_profile_likelihood-1.png" width="100%" />

We only ran the above with one seed. We can draw the curve with using
different seeds to check if this does not change the estimates. We will
likely need to use more samples if the result depends on the seed.

``` r
# assign the seeds we will use along with the grid of seeds and scale parameters
seeds <- 1:4
seeds_n_sigs <- expand.grid(seed = seeds, sig = sigs)

# compute the profile likelihood
ll_terms <- get_pedigree_ll_terms(dat, max_threads = 4L)
pl_curve_res <- Map(function(sig, seed){
  # set the parameters to pass
  beta <- start$beta_no_rng
  sig_sq_log <- 2 * log(sig)
  beta_scaled <- beta * sqrt(1 + sig^2)
  
  # optimize like before but using the fix argument
  opt_out_quick <- pedmod_opt(
    ptr = ll_terms, par = c(beta_scaled, sig_sq_log), maxvls = 1000L, 
    abs_eps = 0, rel_eps = 1e-2, minvls = 100L, use_aprx = TRUE, n_threads = 4L, 
    fix = length(beta) + 1L, seed = seed)
  opt_out <- pedmod_opt(
    ptr = ll_terms, par = c(opt_out_quick$par, sig_sq_log), abs_eps = 0, 
    use_aprx = TRUE, n_threads = 4L, fix = length(beta) + 1L, seed = seed,
    # we changed the parameters
    maxvls = 25000L, rel_eps = 1e-3, minvls = 5000L)
  
  # report to console and return
  message(sprintf("\nLog-likelihood %.5f (%.5f). Estimated parameters:", 
                  -opt_out$value, -opt_out_quick$value))
  message(paste0(capture.output(print(
    c(opt_out$par, Scale = sig))), collapse = "\n"))
  
  list(opt_out_quick = opt_out_quick, opt_out = opt_out)
}, sig = seeds_n_sigs$sig, seed = seeds_n_sigs$seed)
```

The different profile likelihood curves are drawn below:

``` r
# plot the profile likelihood values
pls <- -sapply(pl_curve_res, function(x) x$opt_out$value)
pls <- matrix(pls, length(seeds), dimnames = list(seeds, sigs))
matplot(sigs, t(pls), pch = seq_along(seeds) + 15L, bty = "l", col = "black", 
        xlab = expression(sigma), ylab = "Profile likelihood")
grid()
```

<img src="man/figures/README-draw_mult_pf_curves-1.png" width="100%" />

``` r
# plot the difference from the mean
matplot(sigs, t(pls) - colMeans(pls), pch = seq_along(seeds) + 15L, bty = "l", 
        col = "black", xlab = expression(sigma), 
        ylab = "Diff. profile likelihood")
grid()
```

<img src="man/figures/README-draw_mult_pf_curves-2.png" width="100%" />

We make a small simulation study below where we are interested in the
estimation time and bias.

``` r
# the seeds we will use
seeds <- c(36451989L, 18774630L, 76585289L, 31898455L, 55733878L, 99681114L, 37725150L, 99188448L, 66989159L, 20673587L, 47985954L, 42571905L, 53089211L, 18457743L, 96049437L, 70222325L, 86393368L, 45380572L, 81116968L, 48291155L, 89755299L, 69891073L, 1846862L, 15263013L, 37537710L, 
           25194071L, 14471551L, 38278606L, 55596031L, 5436537L, 75008107L, 83382936L, 50689482L, 71708788L, 52258337L, 23423931L, 61069524L, 24452554L, 32406673L, 14900280L, 24818537L, 59733700L, 82407492L, 95500692L, 62528680L, 88728797L, 9891891L, 36354594L, 69630736L, 41755287L)

# run the simulation study
sim_study <- lapply(seeds, function(s){
  set.seed(s)
  
  # only run the result if it has not been computed
  f <- file.path("cache", "sim_study_simple", paste0("simple-", s, ".RDS"))
  if(!file.exists(f)){
    # simulate the data
    dat <- sim_dat(n_fams = 400L)
    
    # get the starting values
    library(pedmod)
    ll_terms <- get_pedigree_ll_terms(dat, max_threads = 4L)
    ti_start <- system.time(start <- pedmod_start(
      ptr = ll_terms, data = dat, n_threads = 4L))
    start$time <- ti_start
    
    ti_fit <- system.time(
      opt_out <- pedmod_opt(
        ptr = ll_terms, par = start$par, abs_eps = 0, use_aprx = TRUE, 
        n_threads = 4L, 
        maxvls = 25000L, rel_eps = 1e-3, minvls = 5000L))
    opt_out$time <- ti_fit
    
    out <- list(start = start, opt_out = opt_out, 
                ll_no_rng = start$logLik_no_rng)
    saveRDS(out, f)
  }
  
  # report to console and return 
  out <- readRDS(f)
  message(paste0(capture.output(out$opt_out$par), collapse = "\n"))
  
  out
})

# gather the estimates
beta_est <- sapply(sim_study, function(x) head(x$opt_out$par, 3))
sigma_est <- sapply(sim_study, function(x) exp(tail(x$opt_out$par, 1) / 2))

# compute the errors
tmp <- sim_dat(2L)
err <- rbind(beta_est, sigma = sigma_est) - 
  c(attr(tmp, "beta"), sqrt(attr(tmp, "sig_sq")))

# get the bias estimates and the standard errors
rbind(Bias = rowMeans(err), 
      SE   = apply(err, 1, sd) / sqrt(NCOL(err)))
#>      (Intercept) Continuous  Binary   sigma
#> Bias    -0.06431    0.02768 0.03635 0.05534
#> SE       0.04956    0.01684 0.03289 0.03803

# make a box plot
par(mar = c(5, 5, 1, 1))
boxplot(t(err), ylab = "Error")
abline(h = 0, lty = 2)
grid()
```

<img src="man/figures/README-sim_study_simple-1.png" width="100%" />

``` r
# get the average computation times
time_vals <- sapply(sim_study, function(x) {
  keep <- c("opt_out", "start")
  out <- setNames(sapply(x[keep], function(z) z$time["elapsed"]), keep)
  c(out, total = sum(out))
})
summary(t(time_vals))
#>     opt_out         start           total     
#>  Min.   :16.0   Min.   : 1.15   Min.   :17.8  
#>  1st Qu.:20.6   1st Qu.: 1.77   1st Qu.:22.6  
#>  Median :22.6   Median : 2.23   Median :26.1  
#>  Mean   :23.1   Mean   : 3.21   Mean   :26.3  
#>  3rd Qu.:25.8   3rd Qu.: 4.00   3rd Qu.:28.9  
#>  Max.   :32.8   Max.   :10.28   Max.   :36.3
```

### Adding Environmental Effects

As an extension, we can add an environmental effect. The new scale
matrix, the
![C\_{i2}](https://render.githubusercontent.com/render/math?math=C_%7Bi2%7D
"C_{i2}")’s, can be written as:

``` r
C_env <- matrix(0., NROW(fam), NROW(fam))
C_env[c(1:3, 5)   , c(1:3, 5)   ] <- 1
C_env[c(3:4, 7:8) , c(3:4, 7:8) ] <- 1
C_env[c(5:6, 9:10), c(5:6, 9:10)] <- 1

Matrix::Matrix(C_env, sparse = TRUE)
#> 10 x 10 sparse Matrix of class "dsCMatrix"
#>                          
#>  [1,] 1 1 1 . 1 . . . . .
#>  [2,] 1 1 1 . 1 . . . . .
#>  [3,] 1 1 1 1 1 . 1 1 . .
#>  [4,] . . 1 1 . . 1 1 . .
#>  [5,] 1 1 1 . 1 1 . . 1 1
#>  [6,] . . . . 1 1 . . 1 1
#>  [7,] . . 1 1 . . 1 1 . .
#>  [8,] . . 1 1 . . 1 1 . .
#>  [9,] . . . . 1 1 . . 1 1
#> [10,] . . . . 1 1 . . 1 1
```

We assign the new simulation function below but this time we include
only binary covariates:

``` r
# simulates a data set. 
# 
# Args:
#   n_fams: number of families.
#   beta: the fixed effect coefficients.
#   sig_sq: the scale parameters.
sim_dat <- function(n_fams, beta = c(-3, 4), sig_sq = c(2, 1)){
  # setup before the simulations
  Cmat <- 2 * kinship(ped)
  n_obs <- NROW(fam)
  Sig <- diag(n_obs) + sig_sq[1] * Cmat + sig_sq[2] * C_env
  Sig_chol <- chol(Sig)
  
  # simulate the data
  out <- replicate(
    n_fams, {
      # simulate covariates
      X <- cbind(`(Intercept)` = 1, Binary = runif(n_obs) > .9)
      
      # assign the linear predictor + noise
      eta <- drop(X %*% beta) + drop(rnorm(n_obs) %*% Sig_chol)
      
      # return the list in the format needed for the package
      list(y = as.numeric(eta > 0), X = X, scale_mats = list(
        Genetic = Cmat, Environment = C_env))
    }, simplify = FALSE)
  
  # add attributes with the true values and return 
  attributes(out) <- list(beta = beta, sig_sq = sig_sq)
  out
}
```

The model is

  
![\\begin{align\*}&#10; Y\_{ij} &= \\begin{cases} 1 & \\beta\_0 +
\\beta\_1 B\_{ij} + E\_{ij} + G\_{ij} + R\_{ij} \> 0 \\\\ 0 &
\\text{otherwise} \\end{cases} \\\\&#10; X\_{ij} &\\sim N(0, 1)
\\\\&#10; B\_{ij} &\\sim \\text{Bin}(0.5, 1) \\\\&#10; (G\_{i1}, \\dots
G\_{in\_{i}})^\\top &\\sim N^{(n\_i)}(\\vec 0, \\sigma^2\_G C\_{i1})
\\\\&#10;(E\_{i1}, \\dots E\_{in\_{i}})^\\top &\\sim N^{(n\_i)}(\\vec 0,
\\sigma^2\_E C\_{i2}) \\\\&#10; R\_{ij} &\\sim
N(0, 1)\\end{align\*}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign%2A%7D%0A%20Y_%7Bij%7D%20%26%3D%20%5Cbegin%7Bcases%7D%201%20%26%20%5Cbeta_0%20%2B%20%5Cbeta_1%20B_%7Bij%7D%20%2B%20E_%7Bij%7D%20%2B%20G_%7Bij%7D%20%2B%20R_%7Bij%7D%20%3E%200%20%5C%5C%200%20%26%20%5Ctext%7Botherwise%7D%20%5Cend%7Bcases%7D%20%5C%5C%0A%20X_%7Bij%7D%20%26%5Csim%20N%280%2C%201%29%20%5C%5C%0A%20B_%7Bij%7D%20%26%5Csim%20%5Ctext%7BBin%7D%280.5%2C%201%29%20%5C%5C%0A%20%28G_%7Bi1%7D%2C%20%5Cdots%20G_%7Bin_%7Bi%7D%7D%29%5E%5Ctop%20%26%5Csim%20N%5E%7B%28n_i%29%7D%28%5Cvec%200%2C%20%5Csigma%5E2_G%20C_%7Bi1%7D%29%20%5C%5C%0A%28E_%7Bi1%7D%2C%20%5Cdots%20E_%7Bin_%7Bi%7D%7D%29%5E%5Ctop%20%26%5Csim%20N%5E%7B%28n_i%29%7D%28%5Cvec%200%2C%20%5Csigma%5E2_E%20C_%7Bi2%7D%29%20%5C%5C%0A%20R_%7Bij%7D%20%26%5Csim%20N%280%2C%201%29%5Cend%7Balign%2A%7D
"\\begin{align*}
 Y_{ij} &= \\begin{cases} 1 & \\beta_0 + \\beta_1 B_{ij} + E_{ij} + G_{ij} + R_{ij} \> 0 \\\\ 0 & \\text{otherwise} \\end{cases} \\\\
 X_{ij} &\\sim N(0, 1) \\\\
 B_{ij} &\\sim \\text{Bin}(0.5, 1) \\\\
 (G_{i1}, \\dots G_{in_{i}})^\\top &\\sim N^{(n_i)}(\\vec 0, \\sigma^2_G C_{i1}) \\\\
(E_{i1}, \\dots E_{in_{i}})^\\top &\\sim N^{(n_i)}(\\vec 0, \\sigma^2_E C_{i2}) \\\\
 R_{ij} &\\sim N(0, 1)\\end{align*}")  

where
![C\_{i1}](https://render.githubusercontent.com/render/math?math=C_%7Bi1%7D
"C_{i1}") is two times the kinship matrix,
![C\_{i2}](https://render.githubusercontent.com/render/math?math=C_%7Bi2%7D
"C_{i2}") is singular matrix for the environment effect, and
![B\_{ij}](https://render.githubusercontent.com/render/math?math=B_%7Bij%7D
"B_{ij}") is an observed covariate. In this case, we exploit that some
of log marginal likelihood terms are identical. That is, some of the
combinations of pedigrees, covariates, and outcomes match. Therefor, we
can use the `cluster_weights` arguments to reduce the computation time:

``` r
# simulate a data set
set.seed(27107390)
dat <- sim_dat(n_fams = 1000L)

# compute the log marginal likelihood by not using that some of the log marginal 
# likelihood terms are identical
beta_true   <- attr(dat, "beta")
sig_sq_true <- attr(dat, "sig_sq")

library(pedmod)
ll_terms <- get_pedigree_ll_terms(dat, max_threads = 4L)
system.time(ll_res <- eval_pedigree_ll(
  ll_terms, c(beta_true, log(sig_sq_true)), maxvls = 100000L, abs_eps = 0, 
  rel_eps = 1e-3, minvls = 2500L, use_aprx = TRUE, n_threads = 4))
#>    user  system elapsed 
#>   2.108   0.001   0.579
system.time(grad_res <- eval_pedigree_grad(
  ll_terms, c(beta_true, log(sig_sq_true)), maxvls = 100000L, abs_eps = 0, 
  rel_eps = 1e-3, minvls = 2500L, use_aprx = TRUE, n_threads = 4))
#>    user  system elapsed 
#>  65.471   0.005  16.427

# find the duplicated combinations of pedigrees, covariates, and outcomes. One 
# likely needs to change this code if the pedigrees are not identical but are 
# identical if they are permuted. In this case, the code below will miss 
# identical log marginal likelihood terms
dat_unqiue <- dat[!duplicated(dat)]
attributes(dat_unqiue) <- attributes(dat)
length(dat_unqiue) # number of unique terms
#> [1] 416

# get the weights. This can be written in a much more efficient way
c_weights <- sapply(dat_unqiue, function(x)
  sum(sapply(dat, identical, y = x)))

# get the C++ object and show that the computation time is reduced
ll_terms <- get_pedigree_ll_terms(dat_unqiue, max_threads = 4L)

system.time(ll_res_fast <- eval_pedigree_ll(
  ll_terms, c(beta_true, log(sig_sq_true)), maxvls = 100000L, abs_eps = 0, 
  rel_eps = 1e-3, minvls = 2500L, use_aprx = TRUE, n_threads = 4, 
  cluster_weights = c_weights))
#>    user  system elapsed 
#>   1.169   0.000   0.302
system.time(grad_res_fast <- eval_pedigree_grad(
  ll_terms, c(beta_true, log(sig_sq_true)), maxvls = 100000L, abs_eps = 0, 
  rel_eps = 1e-3, minvls = 2500L, use_aprx = TRUE, n_threads = 4, 
  cluster_weights = c_weights))
#>    user  system elapsed 
#>  28.323   0.000   7.105

# show that we get the same (up to a Monte Carlo error)
print(c(redundant = ll_res, fast = ll_res_fast), digits = 6)
#> redundant      fast 
#>  -2632.06  -2632.14
rbind(redundant = grad_res, fast = grad_res_fast)
#>             [,1]   [,2]  [,3]  [,4]
#> redundant -1.854 -5.201 3.434 1.079
#> fast      -1.852 -5.198 3.384 1.069
rm(dat) # will not need this anymore

# find the starting values
start <- pedmod_start(ptr = ll_terms, data = dat_unqiue, 
                      cluster_weights = c_weights)

# optimize
system.time(
  opt_out_quick <- pedmod_opt(
    ptr = ll_terms, par = start$par, abs_eps = 0, use_aprx = TRUE, 
    n_threads = 4L,  cluster_weights = c_weights,
    maxvls = 5000L, rel_eps = 1e-2, minvls = 500L))
#>    user  system elapsed 
#>  15.794   0.003   3.949
system.time(
  opt_out <- pedmod_opt(
    ptr = ll_terms, par = opt_out_quick$par, abs_eps = 0, use_aprx = TRUE, 
    n_threads = 4L,  cluster_weights = c_weights,
    # we changed the parameters
    maxvls = 25000L, rel_eps = 1e-3, minvls = 5000L))
#>    user  system elapsed 
#>  104.49    0.00   26.84
```

The results are shown below:

``` r
# parameter estimates versus the truth
rbind(opt_out       = head(opt_out$par, -2), 
      opt_out_quick = head(start  $par, -2), 
      truth         = attr(dat_unqiue, "beta"))
#>               (Intercept) Binary
#> opt_out            -2.931  3.876
#> opt_out_quick      -2.777  3.669
#> truth              -3.000  4.000
rbind(opt_out       = exp(tail(opt_out$par, 2)), 
      opt_out_quick = exp(tail(start  $par, 2)), 
      truth         = attr(dat_unqiue, "sig_sq"))
#>                           
#> opt_out       1.856 0.9667
#> opt_out_quick 1.555 0.8821
#> truth         2.000 1.0000

# log marginal likelihoods
print( start  $logLik_est, digits = 8)  # this is unreliably/imprecise
#> [1] -2632.1749
print(-opt_out$value     , digits = 8)
#> [1] -2631.9478
```

We can make a 2D profile likelihood curve as follows:

``` r
# get the values at which we evaluate the profile likelihood
rg <- Map(function(est, truth)
  range(exp(est / 2) * c(.8, 1.25), truth), 
  est = tail(opt_out$par, 2), truth = sqrt(attr(dat_unqiue, "sig_sq")))

sig_vals1 <- seq(rg[[1]][1], rg[[1]][2], length.out = 5)
sig_vals2 <- seq(rg[[2]][1], rg[[2]][2], length.out = 5)
sigs <- expand.grid(sigma1 = sig_vals1,
                    sigma2 = sig_vals2)

# function to compute the profile likelihood. 
# 
# Args:
#   fix: indices óf parameters to fix. 
#   fix_val: values of the fixed parameters.
#   sig_start: starting values for the scale parameters.
ll_terms <- get_pedigree_ll_terms(dat_unqiue, max_threads = 4L)
pl_curve_func <- function(fix, fix_val, 
                          sig_start = exp(tail(opt_out$par, 2) / 2)){
  # get the fixed indices of the fixed parameters
  beta = start$beta_no_rng
  is_fix_beta <- fix <= length(beta)
  fix_beta <- fix[is_fix_beta]
  is_fix_sigs <- fix >  length(beta)
  fix_sigs <- fix[is_fix_sigs]
  
  # set the parameters to pass
  sig <- sig_start
  if(length(fix_sigs) > 0)
    sig[fix_sigs - length(beta)] <- fix_val[is_fix_sigs]
  
  # re-scale beta and setup the sigma argument to pass
  sig_sq_log <- 2 * log(sig)
  beta_scaled <- beta * sqrt(1 + sum(sig^2))
  
  # setup the parameter vector
  fix_par <- c(beta_scaled, sig_sq_log)
  if(length(fix_beta) > 0)
    fix_par[fix_beta] <- fix_val[is_fix_beta]
  
  # optimize like before but using the fix argument
  opt_out_quick <- pedmod_opt(
    ptr = ll_terms, par = fix_par, maxvls = 5000L, abs_eps = 0, 
    rel_eps = 1e-2, minvls = 500L, use_aprx = TRUE, n_threads = 4L, 
    fix = fix, cluster_weights = c_weights)
  
  # notice that pedmod_opt only returns a subset of the parameters. These are 
  # the parameters that have been optimized over
  par_new <- fix_par
  par_new[-fix] <- opt_out_quick$par
  opt_out <- pedmod_opt(
    ptr = ll_terms, par = par_new, abs_eps = 0, 
    use_aprx = TRUE, n_threads = 4L, fix = fix,
    cluster_weights = c_weights,
    # we changed the parameters
    maxvls = 25000L, rel_eps = 1e-3, minvls = 5000L)
  
  # report to console and return
  message(sprintf("\nLog-likelihood %.5f (%.5f). Estimated parameters:", 
                  -opt_out$value, -opt_out_quick$value))
  message(paste0(capture.output(print(
    c(`non-fixed` = opt_out$par, fixed = fix_par[fix]))), collapse = "\n"))
  
  list(opt_out_quick = opt_out_quick, opt_out = opt_out)
}

# compute the profile likelihood
pl_curve_res <- Map(
  function(sig1, sig2) pl_curve_func(fix = 0:1 + length(opt_out$par) - 1L, 
                                     fix_val = c(sig1, sig2)), 
  sig1 = sigs$sigma1, sig2 = sigs$sigma2)
```

``` r
par(mfcol = c(2, 2), mar = c(1, 1, 1, 1))
pls <- -sapply(pl_curve_res, function(x) x$opt_out$value)
for(i in 1:3 - 1L)
  persp(sig_vals1, sig_vals2, matrix(pls, length(sig_vals1)), 
        xlab = "\nGenetic", ylab = "\nEnvironment", 
        zlab = "\n\nProfile likelihood", theta = 65 + i * 90, 
        ticktype = "detailed")
```

<img src="man/figures/README-draw_simple_w_ev_ex_profile_likelihood-1.png" width="100%" />

We may just be interested creating two profile likelihood curves for
each of the scale parameters. This can be done as follows:

``` r
# First we compute data for the two profile likelihood curves staring with the
# curve for the additive genetic effect
sigs_genetic <- exp(opt_out$par[3] / 2)
sigs_genetic <- c(sigs_genetic,
                  seq(sigs_genetic - .5, sigs_genetic + .7, length.out = 15))
sigs_genetic <- sort(sigs_genetic)

pl_genetic <- lapply(sigs_genetic, pl_curve_func, fix = 
                       length(start$beta_no_rng) + 1L)

# then we compute the curve for the environmental effect
sigs_env <- exp(opt_out$par[4] / 2)
sigs_env <- c(sigs_env, seq(sigs_env - .3, sigs_env + .5, length.out = 15))
sigs_env <- sort(sigs_env)

pl_env <- lapply(sigs_env, pl_curve_func, fix = 
                       length(start$beta_no_rng) + 2L)
```

``` r
# we create the plots starting with the additive genetic effect  
par(mar = c(5, 5, 1, 1))
pls <- -sapply(pl_genetic, function(x) x$opt_out$value)
plot(sigs_genetic, pls, bty = "l",
     pch = 16, xlab = bquote(paste(sigma['G'])), 
     ylab = "Profile likelihood")
grid()
smooth_est <- smooth.spline(sigs_genetic, pls)
lines(predict(smooth_est, seq(min(sigs_genetic), max(sigs_genetic), 
                              length.out = 100)))
abline(v = exp(opt_out$par[3] / 2), lty = 2) # the estimate
abline(v = sqrt(attr(dat_unqiue, "sig_sq")[1]), lty = 3) # the true value
alpha <- .05
crit_val <- qchisq(1 - alpha, 1)
abline(h = -opt_out$value - crit_val / 2, lty = 3) # mark the critical value
```

<img src="man/figures/README-plot_env_pl_curves-1.png" width="100%" />

``` r
# then we create the plot for the environmental effect 
par(mar = c(5, 5, 1, 1))
pls <- -sapply(pl_env, function(x) x$opt_out$value)
plot(sigs_env, pls, bty = "l",
     pch = 16, xlab = bquote(paste(sigma['E'])), 
     ylab = "Profile likelihood")
grid()
smooth_est <- smooth.spline(sigs_env, pls)
lines(predict(smooth_est, seq(min(sigs_env), max(sigs_env), length.out = 100)))
abline(v = exp(opt_out$par[4] / 2), lty = 2) # the estimate
abline(v = sqrt(attr(dat_unqiue, "sig_sq")[2]), lty = 3) # the true value
abline(h = -opt_out$value - crit_val / 2, lty = 3) # mark the critical value
```

<img src="man/figures/README-plot_env_pl_curves-2.png" width="100%" />

We make a small simulation study below where we are interested in the
estimation time and bias.

``` r
# the seeds we will use
seeds <- c(36451989L, 18774630L, 76585289L, 31898455L, 55733878L, 99681114L, 37725150L, 99188448L, 66989159L, 20673587L, 47985954L, 42571905L, 53089211L, 18457743L, 96049437L, 70222325L, 86393368L, 45380572L, 81116968L, 48291155L, 89755299L, 69891073L, 1846862L, 15263013L, 37537710L, 
           25194071L, 14471551L, 38278606L, 55596031L, 5436537L, 75008107L, 83382936L, 50689482L, 71708788L, 52258337L, 23423931L, 61069524L, 24452554L, 32406673L, 14900280L, 24818537L, 59733700L, 82407492L, 95500692L, 62528680L, 88728797L, 9891891L, 36354594L, 69630736L, 41755287L)

# run the simulation study
sim_study <- lapply(seeds, function(s){
  set.seed(s)
  
  # only run the result if it has not been computed
  f <- file.path("cache", "sim_study_simple_w_env", 
                 paste0("simple-", s, ".RDS"))
  if(!file.exists(f)){
    # simulate the data
    dat <- sim_dat(n_fams = 1000L)
    
    # get the weighted data set
    dat_unqiue <- dat[!duplicated(dat)]
    attributes(dat_unqiue) <- attributes(dat)
    c_weights <- sapply(dat_unqiue, function(x)
      sum(sapply(dat, identical, y = x)))
    rm(dat)
    
    # get the starting values
    library(pedmod)
    ll_terms <- get_pedigree_ll_terms(dat_unqiue, max_threads = 4L)
    ti_start <- system.time(start <- pedmod_start(
      ptr = ll_terms, data = dat_unqiue, n_threads = 4L, 
      cluster_weights = c_weights))
    start$time <- ti_start
    
    # fit the model
    ti_quick <- system.time(
      opt_out_quick <- pedmod_opt(
        ptr = ll_terms, par = start$par, maxvls = 5000L, abs_eps = 0, 
        rel_eps = 1e-2, minvls = 500L, use_aprx = TRUE, n_threads = 4L, 
        cluster_weights = c_weights))
    opt_out_quick$time <- ti_quick
    
    ti_slow <- system.time(
      opt_out <- pedmod_opt(
        ptr = ll_terms, par = opt_out_quick$par, abs_eps = 0, use_aprx = TRUE, 
        n_threads = 4L, cluster_weights = c_weights,
        # we changed the parameters
        maxvls = 25000L, rel_eps = 1e-3, minvls = 5000L))
    opt_out$time <- ti_slow
    
    out <- list(start = start, opt_out = opt_out, opt_out_quick = opt_out_quick, 
                opt_out = opt_out, ll_no_rng = start$logLik_no_rng)
    saveRDS(out, f)
  }
  
  # report to console and return 
  out <- readRDS(f)
  message(paste0(capture.output(out$opt_out$par), collapse = "\n"))
  
  out
})

# gather the estimates
tmp <- sim_dat(2L)
beta_est <- sapply(sim_study, function(x) head(x$opt_out$par, 2))
sigma_est <- sapply(sim_study, function(x) exp(tail(x$opt_out$par, 2) / 2))
rownames(sigma_est) <- names(tmp[[1L]]$scale_mats)

# compute the errors
err <- rbind(beta_est, sigma_est) - 
  c(attr(tmp, "beta"), sqrt(attr(tmp, "sig_sq")))

# get the bias estimates and the standard errors
rbind(Bias = rowMeans(err), 
      SE   = apply(err, 1, sd) / sqrt(NCOL(err)))
#>      (Intercept)  Binary  Genetic Environment
#> Bias    -0.01720 0.02396 0.008509    -0.01165
#> SE       0.04747 0.06488 0.038020     0.01312

# make a box plot
par(mar = c(5, 5, 1, 1))
boxplot(t(err), ylab = "Error")
abline(h = 0, lty = 2)
grid()
```

<img src="man/figures/README-sim_study_simple_w-1.png" width="100%" />

``` r
# get the average computation times
time_vals <- sapply(sim_study, function(x) {
  keep <- c("opt_out", "opt_out_quick", "start")
  out <- setNames(sapply(x[keep], function(z) z$time["elapsed"]), keep)
  c(out, total = sum(out))
})
summary(t(time_vals))
#>     opt_out      opt_out_quick       start          total     
#>  Min.   : 5.08   Min.   : 2.70   Min.   :1.51   Min.   :14.7  
#>  1st Qu.:21.58   1st Qu.: 6.13   1st Qu.:2.06   1st Qu.:33.2  
#>  Median :26.34   Median : 7.29   Median :2.51   Median :40.4  
#>  Mean   :30.36   Mean   : 8.90   Mean   :3.46   Mean   :42.7  
#>  3rd Qu.:34.68   3rd Qu.:12.03   3rd Qu.:4.13   3rd Qu.:45.2  
#>  Max.   :81.74   Max.   :18.93   Max.   :9.74   Max.   :92.8
```

### More Complicated Example

We consider a more complicated example in this section and use some of
the lower level functions in the package as an example. We start by
sourcing a file to get a function to simulate a data set with a maternal
effect and a genetic effect like in Mahjani et al. (2020):

``` r
# source the file to get the simulation function
source(system.file("gen-pedigree-data.R", package = "pedmod"))

# simulate a data set
set.seed(68167102)
dat <- sim_pedigree_data(n_families = 1000L)

# distribution of family sizes
par(mar = c(5, 4, 1, 1))
plot(table(sapply(dat$sim_data, function(x) length(x$y))), 
     xlab = "Family size", ylab = "Number of families", bty = "l")
```

<img src="man/figures/README-source_sim_file-1.png" width="100%" />

``` r
# total number of observations
sum(sapply(dat$sim_data, function(x) length(x$y)))
#> [1] 49734

# average event rate
mean(unlist(sapply(dat$sim_data, `[[`, "y")))
#> [1] 0.2386
```

As Mahjani et al. (2020), we have data families linked by three
generations but we only have data for the last generation. We illustrate
this with the first family in the simulated data set:

``` r
# this is the full family 
library(kinship2)
fam1 <- dat$sim_data[[1L]]
plot(fam1$pedAll)
```

<img src="man/figures/README-one_family-1.png" width="100%" />

``` r
# here is the C matrix for the genetic effect
rev_img <- function(x, ...)
  image(x[, NROW(x):1], ...)
cl <- colorRampPalette(c("Red", "White", "Blue"))(101)

par(mar = c(2, 2, 1, 1))
rev_img(fam1$rel_mat, xaxt = "n", yaxt = "n", col = cl, 
        zlim = c(-1, 1))
```

<img src="man/figures/README-one_family-2.png" width="100%" />

``` r
# the first part of the matrix is given below
with(fam1, rel_mat[seq_len(min(10, NROW(rel_mat))), 
                   seq_len(min(10, NROW(rel_mat)))])
#>        9    10    15    16    17    21    22    28    29    36
#> 9  1.000 0.500 0.125 0.125 0.125 0.000 0.000 0.125 0.125 0.000
#> 10 0.500 1.000 0.125 0.125 0.125 0.000 0.000 0.125 0.125 0.000
#> 15 0.125 0.125 1.000 0.500 0.500 0.125 0.125 0.000 0.000 0.000
#> 16 0.125 0.125 0.500 1.000 0.500 0.125 0.125 0.000 0.000 0.000
#> 17 0.125 0.125 0.500 0.500 1.000 0.125 0.125 0.000 0.000 0.000
#> 21 0.000 0.000 0.125 0.125 0.125 1.000 0.500 0.000 0.000 0.000
#> 22 0.000 0.000 0.125 0.125 0.125 0.500 1.000 0.000 0.000 0.000
#> 28 0.125 0.125 0.000 0.000 0.000 0.000 0.000 1.000 0.500 0.125
#> 29 0.125 0.125 0.000 0.000 0.000 0.000 0.000 0.500 1.000 0.125
#> 36 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.125 0.125 1.000

# here is the C matrix for the maternal effect
rev_img(fam1$met_mat, xaxt = "n", yaxt = "n", col = cl, 
        zlim = c(-1, 1))
```

<img src="man/figures/README-one_family-3.png" width="100%" />

``` r
# each simulated family has such two matrices in addition to a design matrix
# for the fixed effects, X, and a vector with outcomes, y
str(fam1[c("X", "y")])
#> List of 2
#>  $ X: num [1:52, 1:3] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "(Intercept)" "X1" ""
#>  $ y: Named logi [1:52] FALSE TRUE TRUE TRUE FALSE FALSE ...
#>   ..- attr(*, "names")= chr [1:52] "9" "10" "15" "16" ...
```

Then we perform the model estimation:

<!-- knitr::opts_knit$set(output.dir = ".") -->

<!-- knitr::load_cache("est_mod", path = "cache/README-") -->

``` r
# the true parameters are
dat$beta
#> (Intercept)          X1          X2 
#>        -1.0         0.3         0.2
dat$sc # the sigmas squared
#>  Genetic Maternal 
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
#> [1] -26480
(beta <- start_fit$coefficients) # starting values for fixed effects 
#> (Intercept)          X1             
#>     -0.7342      0.2234      0.1349

# start at moderate sized scale parameters
sc <- rep(log(.2), 2)

# check log likelihood at the starting values. First we assign a function 
# to approximate the log likelihood and the gradient
fn <- function(par, seed = 1L, rel_eps = 1e-2, use_aprx = TRUE, 
               n_threads = 4L, indices = NULL, maxvls = 25000L){
  set.seed(seed)
  -eval_pedigree_ll(
    ptr = ll_terms, par = par, maxvls = maxvls, abs_eps = 0, rel_eps = rel_eps, 
    minvls = 1000L, use_aprx = use_aprx, n_threads = n_threads, 
    indices = indices)
}
gr <- function(par, seed = 1L, rel_eps = 1e-2, use_aprx = TRUE, 
               n_threads = 4L, indices = NULL, maxvls = 25000L){
  set.seed(seed)
  out <- -eval_pedigree_grad(
    ptr = ll_terms, par = par, maxvls = maxvls, abs_eps = 0, rel_eps = rel_eps, 
    minvls = 1000L, use_aprx = use_aprx, n_threads = n_threads, 
    indices = indices)
  structure(c(out), value = -attr(out, "logLik"), 
            n_fails = attr(out, "n_fails"))
}

# check output at the starting values
system.time(ll <- -fn(c(beta, sc)))
#>    user  system elapsed 
#>   7.663   0.000   1.954
ll # the log likelihood at the starting values
#> [1] -26042
#> attr(,"n_fails")
#> [1] 0
system.time(gr_val <- gr(c(beta, sc)))
#>    user  system elapsed 
#> 107.690   0.051  27.305
gr_val # the gradient at the starting values
#> [1] 1894.52 -549.92 -235.39   47.16  -47.96
#> attr(,"value")
#> [1] 26042
#> attr(,"n_fails")
#> [1] 721

# standard deviation of the approximation
sd(sapply(1:25, function(seed) fn(c(beta, sc), seed = seed)))
#> [1] 0.09254

# we do the same for the gradient elements but only for a subset of the 
# log marginal likelihood elements
gr_hats <- sapply(1:25, function(seed) gr(c(beta, sc), seed = seed, 
                                          indices = 0:99))
apply(gr_hats, 1, sd)
#> [1] 0.06021 0.09439 0.06420 0.02612 0.02770

# verify the gradient (may not be exactly equal due to MC error)
rbind(numDeriv = numDeriv::grad(fn, c(beta, sc), indices = 0:10), 
      pedmod   = gr(c(beta, sc), indices = 0:10))
#>           [,1]    [,2]  [,3]  [,4]   [,5]
#> numDeriv 28.00 -0.2980 7.415 1.105 -1.071
#> pedmod   27.98 -0.3099 7.395 1.105 -1.076

# optimize the log likelihood approximation
system.time(opt <- optim(c(beta, sc), fn, gr, method = "BFGS"))
#>    user  system elapsed 
#> 4126.00    0.28 1050.03
```

The output from the optimization is shown below:

``` r
print(-opt$value, digits = 8) # the maximum log likelihood
#> [1] -25823.042
opt$convergence               # check convergence
#> [1] 0

# compare the estimated fixed effects with the true values
rbind(truth     = dat$beta, 
      estimated = head(opt$par, length(dat$beta)))
#>           (Intercept)     X1     X2
#> truth          -1.000 0.3000 0.2000
#> estimated      -1.007 0.3075 0.1844

# compare estimated scale parameters with the true values
rbind(truth     = dat$sc, 
      estimated = exp(tail(opt$par, length(dat$sc))))
#>           Genetic Maternal
#> truth      0.5000    0.330
#> estimated  0.5112    0.368
```

### Computation in Parallel

The method scales reasonably well in the number of threads as shown
below:

``` r
library(microbenchmark)
microbenchmark(
  `fn (1 thread)`  = fn(c(beta, sc), n_threads = 1),
  `fn (2 threads)` = fn(c(beta, sc), n_threads = 2),
  `fn (4 threads)` = fn(c(beta, sc), n_threads = 4),
  `gr (1 thread)`  = gr(c(beta, sc), n_threads = 1),
  `gr (2 threads)` = gr(c(beta, sc), n_threads = 2),
  `gr (4 threads)` = gr(c(beta, sc), n_threads = 4),
  times = 1)
#> Unit: seconds
#>            expr     min      lq    mean  median      uq     max neval
#>   fn (1 thread)   7.792   7.792   7.792   7.792   7.792   7.792     1
#>  fn (2 threads)   3.788   3.788   3.788   3.788   3.788   3.788     1
#>  fn (4 threads)   2.086   2.086   2.086   2.086   2.086   2.086     1
#>   gr (1 thread) 101.631 101.631 101.631 101.631 101.631 101.631     1
#>  gr (2 threads)  53.374  53.374  53.374  53.374  53.374  53.374     1
#>  gr (4 threads)  27.365  27.365  27.365  27.365  27.365  27.365     1
```

### Using ADAM

We use stochastic gradient descent with the ADAM method (Kingma and Ba
2015) in this section. We define a function below to apply ADAM and use
it to estimate the model.

<!-- knitr::opts_knit$set(output.dir = ".") -->

<!-- knitr::load_cache("use_adam", path = "cache/README-") -->

``` r
#####
# performs stochastic gradient descent (using ADAM).
#
# Args:
#   par: starting value.
#   gr: function to evaluate the log marginal likelihood.
#   n_clust: number of observation.
#   n_blocks: number of blocks.
#   maxit: maximum number of iteration.
#   seed: seed to use.
#   epsilon, alpha, beta_1, beta_2: ADAM parameters.
#   maxvls: maximum number of samples to draw in each iteration. Thus, it 
#           needs maxit elements.
#   verbose: print output during the estimation.
#   ...: arguments passed to gr.
adam <- function(par, gr, n_clust, n_blocks, maxit = 10L,
                 seed = 1L, epsilon = 1e-8, alpha = .001, beta_1 = .9,
                 beta_2 = .999, maxvls = rep(10000L, maxit), 
                 verbose = FALSE, ...){
  grp_dummy <- (seq_len(n_clust) - 1L) %% n_blocks
  n_par <- length(par)
  m <- v <- numeric(n_par)
  fun_vals <- numeric(maxit)
  estimates <- matrix(NA_real_, n_par, maxit)
  i <- -1L

  for(k in 1:maxit){
    # sample groups
    indices <- sample.int(n_clust, replace = FALSE) - 1L
    blocks <- tapply(indices, grp_dummy, identity, simplify = FALSE)
    
    for(ii in 1:n_blocks){
      i <- i + 1L
      idx_b <- (i %% n_blocks) + 1L
      m_old <- m
      v_old <- v
      res <- gr(par, indices = blocks[[idx_b]], maxvls = maxvls[k])
      fun_vals[(i %/% n_blocks) + 1L] <-
        fun_vals[(i %/% n_blocks) + 1L] + attr(res, "value")
      res <- c(res)

      m <- beta_1 * m_old + (1 - beta_1) * res
      v <- beta_2 * v_old + (1 - beta_2) * res^2

      m_hat <- m / (1 - beta_1^(i + 1))
      v_hat <- v / (1 - beta_2^(i + 1))

      par <- par - alpha * m_hat / (sqrt(v_hat) + epsilon)
    }
    
    if(verbose){
      cat(sprintf("Ended iteration %4d. Running estimate of the function value is: %14.2f\n", 
                  k, fun_vals[k]))
      cat("Parameter estimates are:\n")
      cat(capture.output(print(par)), sep = "\n")
      cat("\n")
    }

    estimates[, k] <- par
  }

  list(par = par, estimates = estimates, fun_vals = fun_vals)
}

#####
# use the function
# assign the maximum number of samples we will use
maxit <- 100L
minvls <- 250L
maxpts <- formals(gr)$maxvls
maxpts_use <- exp(seq(log(2 * minvls), log(maxpts), length.out = maxit))

# show the maximum number of samples we use
par(mar = c(5, 4, 1, 1))
plot(maxpts_use, pch = 16, xlab = "Iteration number", bty = "l",
     ylab = "Maximum number of samples", ylim = range(0, maxpts_use))
```

<img src="man/figures/README-use_adam-1.png" width="100%" />

``` r
set.seed(1)
system.time(
  adam_res <- adam(c(beta, sc), gr = gr, n_clust = length(dat_arg), 
                   n_blocks = 10L, alpha = 1e-2, maxit = maxit, 
                   verbose = FALSE, maxvls = maxpts_use, 
                   minvls = minvls))
#>     user   system  elapsed 
#> 3916.527    0.447 1003.786
```

The result is shown below.

``` r
print(-fn(adam_res$par), digits = 8) # the maximum log likelihood
#> [1] -25823.001
#> attr(,"n_fails")
#> [1] 0

# compare the estimated fixed effects with the true values
rbind(truth             = dat$beta,
      `estimated optim` = head(opt$par     , length(dat$beta)),
      `estimated ADAM`  = head(adam_res$par, length(dat$beta)))
#>                 (Intercept)     X1     X2
#> truth                -1.000 0.3000 0.2000
#> estimated optim      -1.007 0.3075 0.1844
#> estimated ADAM       -1.007 0.3072 0.1860

# compare estimated scale parameters with the true values
rbind(truth             = dat$sc, 
      `estimated optim` = exp(tail(opt$par     , length(dat$sc))), 
      `estimated ADAM`  = exp(tail(adam_res$par, length(dat$sc))))
#>                 Genetic Maternal
#> truth            0.5000   0.3300
#> estimated optim  0.5112   0.3680
#> estimated ADAM   0.5242   0.3635

# could possibly have stopped much earlier maybe. Dashed lines are final 
# estimates
par(mar = c(5, 4, 1, 1))
matplot(t(adam_res$estimates), type = "l", col = "Black", lty = 1, 
        bty = "l", xlab = "Iteration", ylab = "Estimate")
for(s in adam_res$par)
  abline(h = s, lty = 2)
```

<img src="man/figures/README-res_adam-1.png" width="100%" />

### The Multivariate Normal CDF Approximation

We compare the multivariate normal CDF approximation in this section
with the approximation from the mvtnorm package which uses the
implementation by Genz and Bretz (2002). The same algorithm is used but
the version in this package is re-written in C++ and differs slightly.
Moreover, we have implemented an approximation of the standard normal
CDF and its inverse which reduces the computation time as we will show
below.

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
#>        2.830e-05        3.163e-05        3.196e-05
par(mar = c(5, 4, 1, 1))
boxplot(t(sim_res[, "SE", ]))
```

<img src="man/figures/README-show_averge_rel_err-1.png" width="100%" />

The new implementation is faster when the approximation is used:

``` r
rowMeans(sim_res[, "time", ])
#>          mvtnorm mvndst (no aprx) mvndst (w/ aprx) 
#>          0.01833          0.01700          0.01139
par(mar = c(5, 4, 1, 1))
boxplot(t(sim_res[, "time", ]))
```

<img src="man/figures/README-use_new_impl-1.png" width="100%" />

## References

<div id="refs" class="references">

<div id="ref-Genz02">

Genz, Alan, and Frank Bretz. 2002. “Comparison of Methods for the
Computation of Multivariate T Probabilities.” *Journal of Computational
and Graphical Statistics* 11 (4): 950–71.
<https://doi.org/10.1198/106186002394>.

</div>

<div id="ref-Kingma15">

Kingma, Diederik P., and Jimmy Ba. 2015. “Adam: A Method for Stochastic
Optimization.” *CoRR* abs/1412.6980.

</div>

<div id="ref-Liu17">

Liu, Xing-Rong, Yudi Pawitan, and Mark S. Clements. 2017. “Generalized
Survival Models for Correlated Time-to-Event Data.” *Statistics in
Medicine* 36 (29): 4743–62.
<https://doi.org/https://doi.org/10.1002/sim.7451>.

</div>

<div id="ref-Mahjani20">

Mahjani, Behrang, Lambertus Klei, Christina M. Hultman, Henrik Larsson,
Bernie Devlin, Joseph D. Buxbaum, Sven Sandin, and Dorothy E. Grice.
2020. “Maternal Effects as Causes of Risk for Obsessive-Compulsive
Disorder.” *Biological Psychiatry* 87 (12): 1045–51.
<https://doi.org/https://doi.org/10.1016/j.biopsych.2020.01.006>.

</div>

<div id="ref-Pawitan04">

Pawitan, Y., M. Reilly, E. Nilsson, S. Cnattingius, and P. Lichtenstein.
2004. “Estimation of Genetic and Environmental Factors for Binary Traits
Using Family Data.” *Statistics in Medicine* 23 (3): 449–65.
<https://doi.org/https://doi.org/10.1002/sim.1603>.

</div>

</div>
