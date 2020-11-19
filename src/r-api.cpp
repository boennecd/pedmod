#include "pedigree-ll.h"
#include "ped-mem.h"

//' Multivariate Normal Distribution CDF
//' @description
//' Provides an approximation of the multivariate normal distribution CDF
//' over a hyperrectangle.
//'
//' @param lower numeric vector with lower bounds.
//' @param upper numeric vector with upper bounds.
//' @param mu numeric vector with means.
//' @param sigma covariance matrix.
//' @param maxvls maximum number of samples in the approximation.
//' @param abs_eps absolute convergence threshold.
//' @param rel_eps relative convergence threshold.
//' @param minvls minimum number of samples. Negative values provides a
//' default which depends on the dimension of the integration.
//' @param do_reorder \code{TRUE} if a heuristic variable reordering should
//' be used. \code{TRUE} is likely the best value.
//' @param use_aprx \code{TRUE} if a less precise approximation of
//' \code{\link{pnorm}} and \code{\link{qnorm}} should be used. This may
//' reduce the computation time while not affecting the result much.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector mvndst
  (arma::vec const &lower, arma::vec const &upper, arma::vec const &mu,
   arma::mat const &sigma, unsigned const maxvls = 25000,
   double const abs_eps = .001, double const rel_eps = 0,
   int minvls = -1, bool const do_reorder = true,
   bool const use_aprx = false){
  arma::uword const n = lower.n_elem;
  if(upper.n_elem != n)
    throw std::invalid_argument("mvndst: invalid upper");
  if(mu.n_elem != n)
    throw std::invalid_argument("mvndst: invalid mu");
  if(sigma.n_cols != n or sigma.n_rows != n)
    throw std::invalid_argument("mvndst: invalid sigma");
  if(!std::isfinite(abs_eps) or !std::isfinite(rel_eps))
    throw std::invalid_argument("mvndst: invalid abs_eps or rel_eps");

  if(minvls < 0)
    minvls = pedmod::default_minvls(lower.n_elem);

  if(maxvls < static_cast<unsigned>(minvls) or maxvls < 1)
    throw std::invalid_argument("mvndst: invalid maxvls");

  pedmod::likelihood func;
  parallelrng::set_rng_seeds(1);

  pedmod::cdf<pedmod::likelihood>::set_cache(lower.n_elem, 1);
  pedmod::likelihood::set_cache(lower.n_elem, 1);
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

namespace {

struct pedigree_terms {
  unsigned const max_threads;
  std::vector<pedmod::pedigree_ll_term > terms;

  pedigree_terms(Rcpp::List data, unsigned const max_threads):
    max_threads(std::max(1U, max_threads)) {
    terms.reserve(data.size());
    for(auto x : data){
      Rcpp::List xl(static_cast<SEXP>(x)),
      s_mats(static_cast<SEXP>(xl["scale_mats"]));

      arma::mat const X = Rcpp::as<arma::mat>(xl["X"]);
      arma::vec const y = Rcpp::as<arma::vec>(xl["y"]);
      if(y.n_elem > 1000 or y.n_elem < 1)
        throw std::invalid_argument(
            "pedigree_terms: Either dimension zero or dimension greater than 1000");

      std::vector<arma::mat> scale_mats;
      scale_mats.reserve(s_mats.size());
      for(auto &s : s_mats)
        scale_mats.emplace_back(Rcpp::as<arma::mat>(s));

      terms.emplace_back(X, y, scale_mats, max_threads);
    }

    // checks
    if(terms.size() < 1)
      throw std::invalid_argument("pedigree_terms: no terms");
    int const n_fix = terms[0].n_fix_effect,
      n_scales = terms[0].l_factor.scale_mats.size();
    for(auto &tr : terms){
      if(tr.n_fix_effect != n_fix)
        throw std::invalid_argument("pedigree_terms: number of fixed effects do not match");
      if(tr.l_factor.scale_mats.size() != static_cast<size_t>(n_scales))
        throw std::invalid_argument("pedigree_terms: number of scale matrices do not match");
    }
  }
};

Rcpp::IntegerVector get_indices
  (Rcpp::Nullable<Rcpp::IntegerVector> indices,
   pedigree_terms const &terms){
  if(indices.isNotNull())
    return Rcpp::IntegerVector(indices);

  Rcpp::IntegerVector out(terms.terms.size());
  for(int i = 0; i < out.size(); ++i)
    out[i] = i;
  return out;
}

static pedmod::cache_mem<double> r_mem;
} // namespace

//' Get C++ Object for Log Marginal Likelihood Approximations
//'
//' @description
//' Constructs an object needed for \code{\link{eval_pedigree_ll}} and
//' \code{\link{eval_pedigree_grad}}.
//'
//' @param data \code{\link{list}} where each element is a list for a cluster
//' with an:
//' \itemize{
//'   \item{\code{"X"}}{ element with the design matrix for the fixed effect,}
//'   \item{\code{"y"}}{ element with the zero-one outcomes, and}
//'   \item{\code{"scale_mats"}}{ element with a list where each element is a
//' scale/correlation matrix for a particular type of effect.}
//' }
//' @param max_threads maximum number of threads to use.
//'
//' @export
// [[Rcpp::export]]
SEXP get_pedigree_ll_terms(Rcpp::List data, unsigned const max_threads){
  return Rcpp::XPtr<pedigree_terms>(new pedigree_terms(data, max_threads));
}

//' @rdname eval_pedigree
//'
//' @title Approximate the Log Marginal Likelihood
//' @description
//' Approximate the log marginal likelihood and the derivatives with
//' respect to the model parameters.
//'
//' @param ptr object from \code{\link{get_pedigree_ll_terms}}.
//' @param par numeric vector with fixed effect coefficients and log scale
//' parameters.
//' @param indices zero-based vector with indices of which log marginal
//' likelihood terms to include. Use \code{NULL} if all indices should be
//' used.
//' @param n_threads number of threads to use.
//'
//' @inheritParams get_pedigree_ll_terms
//' @inheritParams mvndst
//'
//' @return \code{eval_pedigree_ll}:
//' a scalar with the log marginal likelihood approximation.
//'
//' @export
// [[Rcpp::export]]
double eval_pedigree_ll
  (SEXP ptr, arma::vec par, int const maxvls,
   double const abs_eps, double const rel_eps,
   Rcpp::Nullable<Rcpp::IntegerVector> indices = R_NilValue,
   int const minvls = -1, bool const do_reorder = true,
   bool const use_aprx = false, unsigned n_threads = 1L){
#ifndef _OPENMP
  n_threads = 1L;
#endif

  Rcpp::XPtr<pedigree_terms> terms_ptr(ptr);
  std::vector<pedmod::pedigree_ll_term > &terms = terms_ptr->terms;
  n_threads = std::min(n_threads, terms_ptr->max_threads);

  parallelrng::set_rng_seeds(n_threads);

  // checks
  int const n_fix = terms[0].n_fix_effect,
         n_scales = terms[0].l_factor.scale_mats.size();
  if(static_cast<int>(par.size()) != n_fix + n_scales)
    throw std::invalid_argument("eval_pedigree_ll: invalid par parameter");

  // transform scale parameters
  for(int i = n_fix; i < n_fix + n_scales; ++i)
    par[i] = std::exp(par[i]);

  r_mem.set_n_mem(1, n_threads);

  // compute
  auto all_idx = get_indices(indices, *terms_ptr);
  int const * idx = &all_idx[0];

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
{
#endif
  double *wmem = r_mem.get_mem();
  *wmem = 0;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for(int i = 0; i < all_idx.size(); ++i){
    if(idx[i] >= static_cast<int>(terms.size()))
      continue;
    *wmem += terms.at(idx[i]).fn(
      &par[0], maxvls, abs_eps, rel_eps, minvls, do_reorder, use_aprx);
  }
#ifdef _OPENMP
}
#endif

  double out(0.);
  for(unsigned i = 0; i < n_threads; ++i)
    out += *r_mem.get_mem(i);

  return out;
}

//' @rdname eval_pedigree
//'
//' @return \code{eval_pedigree}: a vector with the derivatives with
//' respect to \code{par}. An attribute called \code{"logLik"} contains the
//' log marginal likelihood approximation.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector eval_pedigree_grad
  (SEXP ptr, arma::vec par, int const maxvls,
   double const abs_eps, double const rel_eps,
   Rcpp::Nullable<Rcpp::IntegerVector> indices = R_NilValue,
   int const minvls = -1, bool const do_reorder = true,
   bool const use_aprx = false, unsigned n_threads = 1L){
#ifndef _OPENMP
  n_threads = 1L;
#endif

  Rcpp::XPtr<pedigree_terms> terms_ptr(ptr);
  std::vector<pedmod::pedigree_ll_term > &terms = terms_ptr->terms;
  n_threads = std::min(n_threads, terms_ptr->max_threads);

  parallelrng::set_rng_seeds(n_threads);

  // checks
  int const n_fix = terms[0].n_fix_effect,
    n_scales = terms[0].l_factor.scale_mats.size();
  if(static_cast<int>(par.size()) != n_fix + n_scales)
    throw std::invalid_argument("eval_pedigree_ll: invalid par parameter");

  // transform scale parameters
  for(unsigned i = n_fix; i < par.size(); ++i)
    par[i] = std::exp(par[i]);

  r_mem.set_n_mem(1 + par.size(), n_threads);

  // compute
  auto all_idx = get_indices(indices, *terms_ptr);
  int const * idx = &all_idx[0];

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
{
#endif
  double *wmem = r_mem.get_mem();
  std::fill(wmem, wmem + 1 + par.size(), 0.);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for(int i = 0; i < all_idx.size(); ++i){
    if(idx[i] >= static_cast<int>(terms.size()))
      continue;
    *wmem += terms.at(idx[i]).gr(
      &par[0], wmem + 1, maxvls, abs_eps, rel_eps, minvls, do_reorder,
      use_aprx);
  }
#ifdef _OPENMP
}
#endif

  // aggregate the result
  Rcpp::NumericVector grad(n_fix + n_scales);
  double ll(0.);
  for(unsigned i = 0; i < n_threads; ++i){
    double *wmem = r_mem.get_mem(i);
    ll += *wmem;
    for(unsigned j = 0; j <  par.size(); ++j)
      grad[j] += wmem[j + 1];
  }

  // account for exp(...)
  for(int i = n_fix; i < n_fix + n_scales; ++i)
    grad[i] *= par[i];
  grad.attr("logLik") = Rcpp::NumericVector::create(ll);

  return grad;
}
