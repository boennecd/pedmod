#include "norm-cdf-approx.h"

// #include "Rcpp.h"
//
// // [[Rcpp::export(rng = false)]]
// Rcpp::NumericVector test_func(Rcpp::NumericVector x){
//   Rcpp::NumericVector out = Rcpp::clone(x);
//   for(auto &o : out)
//     o = pnorm_approx(o);
//
//   return out;
// }

/*** R
# data for the spline
m_splin <- local({
  n_points <- 210L
  eps <- 1e-10
  x <- seq(qnorm(eps), 0, length.out = n_points)
  splinefun(x = x, y = pnorm(x), method = "monoH.FC")
})

ev <- environment(m_splin)
dput(ev$dx[1]) # h
dput(c(mapply(c, ev$x0, ev$m, ev$y0)))

# test the function
test_vals <- exp(seq(
  log(.Machine$double.eps) * 8, log(1 - .Machine$double.eps) / 8,
  length.out = 10000))
test_vals <- qnorm(test_vals)

truth <- pnorm(test_vals)
res <- test_func(test_vals)
max(abs(truth - res) / truth)
range(truth - res)
all(res >= 0)
min(res)

plot(test_vals, (truth - res) / truth, type = "h")
plot(test_vals, abs((truth - res) / truth), type = "h", log = "y")
plot(test_vals,  truth - res, type = "h")

test_vals <- -test_vals

truth <- pnorm(test_vals)
res <- test_func(test_vals)
max(abs(truth - res) / truth)
range(truth - res)
all(res >= 0)
min(res)

plot(test_vals, (truth - res) / truth, type = "h")
plot(test_vals, abs((truth - res) / truth), type = "h", log = "y")
plot(test_vals,  truth - res, type = "h")

bench::mark(pnorm(test_vals), test_func(test_vals), min_time = 2,
            check = FALSE)

x <- rnorm(10000)
bench::mark(pnorm(x), test_func(x), min_time = 2,
            check = FALSE)
*/
