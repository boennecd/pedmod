#include "cdfarpx.h"

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector mvndst
  (arma::vec const &lower, arma::vec const &upper, arma::vec const &mu,
   arma::mat const &sigma, unsigned const maxvls = 25000,
   double const abs_eps = .001, double const rel_eps = 0,
   int minvls = -1, bool const do_reorder = true,
   bool const use_aprx = false){
  if(minvls < 0)
    minvls = pedmod::default_minvls(lower.n_elem);

  pedmod::likelihood func;
  parallelrng::set_rng_seeds(1);
  auto const out = pedmod::cdf<pedmod::likelihood>(
    func, lower, upper, mu, sigma, do_reorder, use_aprx).approximate(
        maxvls, abs_eps, rel_eps, minvls);

  Rcpp::NumericVector res(1);
  res[0] = out.likelihood;
  res.attr("n_it")   = Rcpp::IntegerVector::create(out.minvls);
  res.attr("inform") = Rcpp::IntegerVector::create(out.inform);
  res.attr("abserr") = Rcpp::NumericVector::create(out.abserr);
  return res;
}
