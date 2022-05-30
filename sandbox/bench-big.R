library(microbenchmark)
dat <- readRDS("dat-big.RDS")
par <- attr(dat, "par")

library(pedmod)
ll_terms <- pedigree_ll_terms(dat, max_threads = 4L)

do_eval <- function(n_threads, what, use_aprx){
  f <- switch(what,
              likelihood = eval_pedigree_ll,
              gradient = eval_pedigree_grad,
              hess = eval_pedigree_hess,
              stop())

  f(ptr = ll_terms, par = par, maxvls = 2500, minvls = 2500,
    abs_eps = 0, rel_eps = 0, use_aprx = use_aprx, n_threads = n_threads)
}

set.seed(1)
bench <- microbenchmark(
  `fn 1` = do_eval(1, "likelihood", TRUE),
  `fn 2` = do_eval(2, "likelihood", TRUE),
  `fn 4` = do_eval(4, "likelihood", TRUE),
  `gr 1` = do_eval(1, "gradient", TRUE),
  `gr 2` = do_eval(2, "gradient", TRUE),
  `gr 4` = do_eval(4, "gradient", TRUE),
  `hess 1` = do_eval(1, "hess", TRUE),
  `hess 2` = do_eval(2, "hess", TRUE),
  `hess 4` = do_eval(4, "hess", TRUE),
  times = 10)
summary(bench)[, c("expr", "lq", "median", "uq")]
