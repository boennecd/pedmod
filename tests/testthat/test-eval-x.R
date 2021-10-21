context("testing eval functions")

test_that("A previous problematic case gives the correct result", {
  dat <- list(list(y = 1L, X = matrix(c(1, 1, 34), 1L),
                   scale_mats = list(matrix(1, 1, 1))))
  ptr <- pedigree_ll_terms(dat, 1L)

  par <- c(-3.64537483404081, 0.06815334989362926, 0.0162241770644372,
           -.487601434379729)

  fn_func <- function(par)
    pnorm(sum(dat[[1L]]$X * par[-4]) / sqrt(1 + exp(par[4])), log.p = TRUE)

  # check that we get the right function value
  fn <- eval_pedigree_ll(ptr, par, maxvls = 1000L, abs_eps = 0,
                         minvls = 100, rel_eps = 1e-3, n_threads = 1L)
  expect_equal(c(fn), fn_func(par), log.p = TRUE)

  fn <- eval_pedigree_ll(ptr, par, maxvls = 1000L, abs_eps = 0,
                         minvls = 100, rel_eps = 1e-3, n_threads = 1L,
                         use_aprx = TRUE)
  expect_equal(c(fn), fn_func(par), log.p = TRUE)

  # dput(numDeriv::grad(fn_func, par))
  true_grad <- c(2.13706257290667, 2.13706257261275, 72.6601274806369, 1.23000330372289)
  gr <- eval_pedigree_grad(ptr, par, maxvls = 1000L, abs_eps = 0,
                           minvls = 100, rel_eps = 1e-3, n_threads = 1L)
  expect_equal(true_grad, c(gr), tolerance = .Machine$double.eps^(4/11))
})
