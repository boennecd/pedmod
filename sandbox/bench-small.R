library(microbenchmark)
dat <- readRDS("dat-small.RDS")
par <- c(attr(dat, "beta"), log(attr(dat, "sig_sq")))

library(pedmod)
ll_terms <- get_pedigree_ll_terms(dat, max_threads = 4L)

do_eval <- function(n_threads, is_grad, use_aprx){
  f <- if(is_grad) eval_pedigree_grad else eval_pedigree_ll
  f(ptr = ll_terms, par = par, maxvls = 2500, minvls = 2500,
    abs_eps = 0, rel_eps = 0, use_aprx = use_aprx, n_threads = n_threads)
}

set.seed(1)
bench <- microbenchmark(
  `fn 1` = do_eval(1, FALSE, TRUE),
  `fn 2` = do_eval(2, FALSE, TRUE),
  `fn 4` = do_eval(4, FALSE, TRUE),
  `gr 1` = do_eval(1, TRUE , TRUE),
  `gr 2` = do_eval(2, TRUE , TRUE),
  `gr 4` = do_eval(4, TRUE , TRUE),
  times = 10)
summary(bench)[, c("expr", "lq", "median", "uq")]
