library(kinship2)

# create the family we will use
fam <- data.frame(id = 1:10, sex = rep(1:2, 5L),
                  father = c(NA, NA, 1L, NA, 1L, NA, 3L, 3L, 5L, 5L),
                  mother = c(NA, NA, 2L, NA, 2L, NA, 4L, 4L, 6L, 6L))

ped <- with(fam, pedigree(id = id, dadid = father, momid = mother, sex = sex))

C_env <- diag(1, NROW(fam))
C_env[c(3, 5), c(3, 5)] <- 1
C_env[c(7:8 ), c(7:8 )] <- 1
C_env[c(9:10), c(9:10)] <- 1

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

set.seed(95635332)
saveRDS(sim_dat(n_fams = 100L), "dat-small.RDS")

# simulate the big data set
rm(list = ls())

# source the file to get the simulation function
source(system.file("gen-pedigree-data.R", package = "pedmod"))

# simulate a data set
set.seed(95635332)
dat <- sim_pedigree_data(n_families = 10L, max_members = 300L)

dat <- structure(
  lapply(dat$sim_data, function(x)
    list(y = as.numeric(x$y), X = x$X,
        scale_mats = list(x$rel_mat, x$met_mat))),
  par = c(dat$beta, dat$sc))
saveRDS(dat, "dat-big.RDS")
