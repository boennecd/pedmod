#include "pedigree-ll.h"

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

//' @export
// [[Rcpp::export]]
SEXP get_pedigree_ll_terms(Rcpp::List data){
  Rcpp::XPtr<std::vector<pedmod::pedigree_ll_term> > terms_ptr(
      new std::vector<pedmod::pedigree_ll_term >());
  std::vector<pedmod::pedigree_ll_term > &terms = *terms_ptr;
  terms.reserve(data.size());

  for(auto x : data){
    Rcpp::List xl(static_cast<SEXP>(x)),
           s_mats(static_cast<SEXP>(xl["scale_mats"]));

    arma::mat const X = Rcpp::as<arma::mat>(xl["X"]);
    arma::vec const y = Rcpp::as<arma::vec>(xl["y"]);
    std::vector<arma::mat> scale_mats;
    scale_mats.reserve(s_mats.size());
    for(auto &s : s_mats)
      scale_mats.emplace_back(Rcpp::as<arma::mat>(s));

    terms.emplace_back(X, y, scale_mats);
  }

  // checks
  if(terms.size() < 1)
    throw std::invalid_argument("get_pedigree_ll_terms: no terms");
  int const n_fix = terms[0].n_fix_effect,
         n_scales = terms[0].l_factor.scale_mats.size();
  for(auto &tr : terms){
    if(tr.n_fix_effect != n_fix)
      throw std::invalid_argument("get_pedigree_ll_terms: number of fixed effects do not match");
    if(tr.l_factor.scale_mats.size() != static_cast<size_t>(n_scales))
      throw std::invalid_argument("get_pedigree_ll_terms: number of scale matrices do not match");
  }

  return terms_ptr;
}

//' @export
// [[Rcpp::export]]
double eval_pedigree_ll
  (SEXP ptr, arma::vec par, int const maxvls,
   double const abs_eps, double const rel_eps, int const minvls = -1,
   bool const do_reorder = true, bool const use_aprx = false,
   unsigned const n_threads = 1L){
  parallelrng::set_rng_seeds(n_threads);
  Rcpp::XPtr<std::vector<pedmod::pedigree_ll_term> > terms_ptr(ptr);
  std::vector<pedmod::pedigree_ll_term > &terms = *terms_ptr;

  // checks
  int const n_fix = terms[0].n_fix_effect,
         n_scales = terms[0].l_factor.scale_mats.size();
  if(static_cast<int>(par.size()) != n_fix + n_scales)
    throw std::invalid_argument("eval_pedigree_ll: invalid par parameter");

  // transform scale parameters
  for(int i = n_fix; i < n_fix + n_scales; ++i)
    par[i] = std::exp(par[i]);

  // compute
  double out(0.);
  for(auto &t : terms)
    out += t.fn(&par[0], maxvls, abs_eps, rel_eps, minvls, do_reorder,
                use_aprx);

  return out;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector eval_pedigree_grad
  (SEXP ptr, arma::vec par, int const maxvls,
   double const abs_eps, double const rel_eps, int const minvls = -1,
   bool const do_reorder = true, bool const use_aprx = false,
   unsigned const n_threads = 1L){
  parallelrng::set_rng_seeds(n_threads);
  Rcpp::XPtr<std::vector<pedmod::pedigree_ll_term> > terms_ptr(ptr);
  std::vector<pedmod::pedigree_ll_term > &terms = *terms_ptr;

  // checks
  int const n_fix = terms[0].n_fix_effect,
    n_scales = terms[0].l_factor.scale_mats.size();
  if(static_cast<int>(par.size()) != n_fix + n_scales)
    throw std::invalid_argument("eval_pedigree_ll: invalid par parameter");

  // transform scale parameters
  for(int i = n_fix; i < n_fix + n_scales; ++i)
    par[i] = std::exp(par[i]);

  // compute
  Rcpp::NumericVector grad(n_fix + n_scales);
  double ll(0.);
  for(auto &t : terms)
    ll += t.gr(&par[0], &grad[0], maxvls, abs_eps, rel_eps, minvls,
               do_reorder, use_aprx);

  for(int i = n_fix; i < n_fix + n_scales; ++i)
    grad[i] *= par[i];
  grad.attr("logLik") = Rcpp::NumericVector::create(ll);

  return grad;
}
