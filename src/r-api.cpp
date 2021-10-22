#include "pedigree-ll.h"
#include "ped-mem.h"
#include "openmp-exception_ptr.h"

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
//' @param method integer with the method to use. Zero yields randomized Korobov
//' lattice rules while one yields scrambled Sobol sequences.
//' @param n_sequences number of randomized quasi-Monte Carlo sequences to use.
//' More samples yields a better estimate of the error but a worse
//' approximation. Eight is used in the original Fortran code. If one is
//' used then the error will be set to zero because it cannot be estimated.
//'
//' @return
//' An approximation of the CDF. The \code{"n_it"} attribute shows the number of
//' integrand evaluations, the \code{"inform"} attribute is zero if the
//' requested precision is achieved, and the \code{"abserr"} attribute
//' shows 3.5 times the estimated standard error.
//'
//' @examples
//' # simulate covariance matrix and the upper bound
//' set.seed(1)
//' n <- 10L
//' S <- drop(rWishart(1L, 2 * n, diag(n) / 2 / n))
//' u <- rnorm(n)
//'
//' system.time(pedmod_res <- mvndst(
//'   lower = rep(-Inf, n), upper = u, sigma = S, mu = numeric(n),
//'   maxvls = 1e6, abs_eps = 0, rel_eps = 1e-4, use_aprx = TRUE))
//' pedmod_res
//'
//' # compare with mvtnorm
//' if(require(mvtnorm)){
//'   mvtnorm_time <- system.time(mvtnorm_res <- pmvnorm(
//'     upper = u, sigma = S, algorithm = GenzBretz(
//'       maxpts = 1e6, abseps = 0, releps = 1e-4)))
//'   cat("mvtnorm_res:\n")
//'   print(mvtnorm_res)
//'
//'   cat("mvtnorm_time:\n")
//'   print(mvtnorm_time)
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector mvndst
  (arma::vec const &lower, arma::vec const &upper, arma::vec const &mu,
   arma::mat const &sigma, unsigned const maxvls = 25000,
   double const abs_eps = .001, double const rel_eps = 0,
   int minvls = -1, bool const do_reorder = true,
   bool const use_aprx = false, int const method = 0,
   unsigned const n_sequences = 8){
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

  pedmod::cdf<pedmod::likelihood>::alloc_mem(lower.n_elem, 1);
  pedmod::likelihood::alloc_mem(lower.n_elem, 1, n_sequences);
  auto const out = pedmod::cdf<pedmod::likelihood>(
    func, lower, upper, mu, sigma, do_reorder, use_aprx).approximate(
        maxvls, abs_eps, rel_eps, pedmod::get_cdf_methods(method), minvls,
        n_sequences);

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

  pedigree_terms(Rcpp::List data, unsigned const max_threads,
                 unsigned const n_sequences):
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

      terms.emplace_back(X, y, scale_mats, max_threads, n_sequences);
    }

    // checks
    if(terms.size() < 1)
      throw std::invalid_argument("pedigree_terms: no terms");
    unsigned const n_fix = terms[0].n_fix_effect,
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

inline unsigned eval_get_n_threads(unsigned const n_threads,
                                   pedigree_terms const &terms){
#ifndef _OPENMP
  return 1L;
#endif

  if(n_threads > terms.max_threads){
    Rcpp::Function warning("warning");
    warning("Cannot set the number of threads to ", std::to_string(n_threads),
            ". The object is created with a maximum of ",
            std::to_string(terms.max_threads), " threads.");
  }

  return std::min(n_threads, terms.max_threads);
}
} // namespace

//' Get a C++ Object for Log Marginal Likelihood Approximations
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
//' @param n_sequences number of randomized quasi-Monte Carlo sequences to use.
//' More samples yields a better estimate of the error but a worse
//' approximation. Eight is used in the original Fortran code. If one is
//' used then the error will be set to zero because it cannot be estimated.
//'
//' @details
//' An intercept column is not added to the \code{X} matrices
//' like what \code{\link{lm.fit}} and \code{\link{glm.fit}} do.
//' Thus, it is often important that the user adds an intercept column
//' to these matrices as it is hardly ever justified to not include the
//' intercept (the exceptions being e.g. when splines are used which include
//' the intercept and with certain dummy designs).
//'
//' @examples
//' # three families as an example
//' fam_dat <- list(
//'   list(
//'     y = c(FALSE, TRUE, FALSE, FALSE),
//'     X = structure(c(
//'       1, 1, 1, 1, 1.2922654151273, 0.358134905909256, -0.734963997107464,
//'       0.855235473516044, -1.16189500386223, -0.387298334620742,
//'       0.387298334620742, 1.16189500386223),
//'       .Dim = 4:3, .Dimnames = list( NULL, c("(Intercept)", "X1", ""))),
//'     rel_mat = structure(c(
//'       1, 0.5, 0.5, 0.125, 0.5, 1, 0.5, 0.125, 0.5, 0.5,
//'       1, 0.125, 0.125, 0.125, 0.125, 1), .Dim = c(4L, 4L)),
//'     met_mat = structure(c(1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1),
//'                         .Dim = c(4L, 4L))),
//'   list(
//'     y = c(FALSE, FALSE, FALSE),
//'     X = structure(c(
//'       1, 1, 1, -0.0388728997202442, -0.0913782435233639,
//'       -0.0801619722392612, -1, 0, 1), .Dim = c(3L, 3L)),
//'     rel_mat = structure(c(
//'       1, 0.5, 0.125, 0.5, 1, 0.125, 0.125, 0.125, 1), .Dim = c(3L, 3L)),
//'     met_mat = structure(c(
//'       1, 1, 0, 1, 1, 0, 0, 0, 1), .Dim = c(3L, 3L))),
//'   list(
//'     y = c(TRUE, FALSE),
//'     X = structure(c(
//'       1, 1, 0.305275750370738, -1.49482995913648,  -0.707106781186547,
//'       0.707106781186547),
//'       .Dim = 2:3, .Dimnames = list( NULL, c("(Intercept)", "X1", ""))),
//'     rel_mat = structure(c(1, 0.5,  0.5, 1), .Dim = c(2L, 2L)),
//'     met_mat = structure(c(1, 1, 1, 1), .Dim = c(2L,  2L))))
//'
//' # get the data into the format needed for the package
//' dat_arg <- lapply(fam_dat, function(x){
//'   # we need the following for each family:
//'   #   y: the zero-one outcomes.
//'   #   X: the design matrix for the fixed effects.
//'   #   scale_mats: list with the scale matrices for each type of effect.
//'   list(y = as.numeric(x$y), X = x$X,
//'        scale_mats = list(x$rel_mat, x$met_mat))
//' })
//'
//' # get a pointer to the C++ object
//' ptr <- pedigree_ll_terms(dat_arg, max_threads = 1L)
//'
//' @export
// [[Rcpp::export]]
SEXP pedigree_ll_terms(Rcpp::List data, unsigned const max_threads = 1,
                       unsigned const n_sequences = 8){
  return Rcpp::XPtr<pedigree_terms>(
    new pedigree_terms(data, max_threads, n_sequences));
}

// [[Rcpp::export]]
int get_n_scales(SEXP ptr){
  return Rcpp::XPtr<pedigree_terms>(ptr)->terms[0].l_factor.scale_mats.size();
}

// [[Rcpp::export]]
int get_n_terms(SEXP ptr){
  return Rcpp::XPtr<pedigree_terms>(ptr)->terms.size();
}

// [[Rcpp::export("eval_pedigree_ll_cpp")]]
Rcpp::NumericVector eval_pedigree_ll
  (SEXP ptr, arma::vec par, int const maxvls,
   double const abs_eps, double const rel_eps,
   Rcpp::Nullable<Rcpp::IntegerVector> indices = R_NilValue,
   int const minvls = -1, bool const do_reorder = true,
   bool const use_aprx = false, unsigned n_threads = 1L,
   Rcpp::Nullable<Rcpp::NumericVector> cluster_weights = R_NilValue,
   int const method = 0){
  Rcpp::XPtr<pedigree_terms> terms_ptr(ptr);
  std::vector<pedmod::pedigree_ll_term > &terms = terms_ptr->terms;
  n_threads = eval_get_n_threads(n_threads, *terms_ptr);

  parallelrng::set_rng_seeds(n_threads);

  // checks
  int const n_fix = terms[0].n_fix_effect,
         n_scales = terms[0].l_factor.scale_mats.size();
  if(static_cast<int>(par.size()) != n_fix + n_scales)
    throw std::invalid_argument(
        "eval_pedigree_ll: invalid par argument. Had " +
          std::to_string(par.size()) + " elements but should have " +
          std::to_string(n_fix + n_scales) + ".");

  if(maxvls < minvls or maxvls < 1)
    throw std::invalid_argument("mvndst: invalid maxvls");

  // get potential weights
  arma::vec c_weights;
  if(cluster_weights.isNotNull()){
    Rcpp::NumericVector r_weights(cluster_weights);
    if(static_cast<size_t>(r_weights.size()) != terms.size())
      throw std::invalid_argument(
          "invalid size of cluster_weights. Should have length " +
            std::to_string(terms.size()) +  " had length " +
            std::to_string(r_weights.size()) + ".");
    c_weights = r_weights;
  }
  bool const has_weights = c_weights.size() > 0;

  // transform scale parameters
  for(int i = n_fix; i < n_fix + n_scales; ++i)
    par[i] = std::exp(par[i]);

  pedmod::cache_mem<double> r_mem;
  r_mem.set_n_mem(2, n_threads);

  // compute
  auto all_idx = get_indices(indices, *terms_ptr);
  int const * idx = &all_idx[0];

  int n_fails(0);
  openmp_exception_ptr exception_handler;
  pedmod::cdf_methods const meth = pedmod::get_cdf_methods(method);

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
{
#endif
  double *wmem = r_mem.get_mem();
  wmem[0] = 0;
  wmem[1] = 0;

#ifdef _OPENMP
#pragma omp for schedule(static) reduction(+:n_fails)
#endif
  for(int i = 0; i < all_idx.size(); ++i)
    exception_handler.run([&]() -> void {
      if(idx[i] >= static_cast<int>(terms.size()))
        return;
      bool did_fail(false);
      double const w_i = has_weights ? c_weights[idx[i]] : 1;
      if(std::abs(w_i) < std::numeric_limits<double>::epsilon())
        return;

      auto const res = terms.at(idx[i]).fn(
        &par[0], maxvls, abs_eps, rel_eps, minvls, do_reorder, use_aprx,
        did_fail, meth);

      wmem[0] += w_i * res.log_likelihood;
      wmem[1] += w_i * w_i * res.estimator_var;
      n_fails += did_fail;
    });
#ifdef _OPENMP
}
#endif

  exception_handler.rethrow_if_error();

  double out(0.),
     var_est(0.);
  for(unsigned i = 0; i < n_threads; ++i){
    out     += r_mem.get_mem(i)[0];
    var_est += r_mem.get_mem(i)[1];
  }

  Rcpp::NumericVector v_out = Rcpp::NumericVector::create(out);
  v_out.attr("n_fails") = Rcpp::IntegerVector::create(n_fails);
  v_out.attr("std"    ) = Rcpp::NumericVector::create(std::sqrt(var_est));
  return v_out;
}


// [[Rcpp::export("eval_pedigree_grad_cpp")]]
Rcpp::NumericVector eval_pedigree_grad
  (SEXP ptr, arma::vec par, int const maxvls,
   double const abs_eps, double const rel_eps,
   Rcpp::Nullable<Rcpp::IntegerVector> indices = R_NilValue,
   int const minvls = -1, bool const do_reorder = true,
   bool const use_aprx = false, unsigned n_threads = 1L,
   Rcpp::Nullable<Rcpp::NumericVector> cluster_weights = R_NilValue,
   int const method = 0){
  Rcpp::XPtr<pedigree_terms> terms_ptr(ptr);
  std::vector<pedmod::pedigree_ll_term > &terms = terms_ptr->terms;
  n_threads = eval_get_n_threads(n_threads, *terms_ptr);

  parallelrng::set_rng_seeds(n_threads);

  // checks
  int const n_fix = terms[0].n_fix_effect,
         n_scales = terms[0].l_factor.scale_mats.size();
  if(static_cast<int>(par.size()) != n_fix + n_scales)
    throw std::invalid_argument(
        "eval_pedigree_grad: invalid par argument. Had " +
          std::to_string(par.size()) + " elements but should have " +
          std::to_string(n_fix + n_scales) + ".");

  // get potential weights
  arma::vec c_weights;
  if(cluster_weights.isNotNull()){
    Rcpp::NumericVector r_weights(cluster_weights);
    if(static_cast<size_t>(r_weights.size()) != terms.size())
      throw std::invalid_argument(
          "invalid size of cluster_weights. Should have length " +
            std::to_string(terms.size()) +  " had length " +
            std::to_string(r_weights.size()) + ".");
    c_weights = r_weights;
  }
  bool const has_weights = c_weights.size() > 0;

  // transform scale parameters
  for(unsigned i = n_fix; i < par.size(); ++i)
    par[i] = std::exp(par[i]);

  pedmod::cache_mem<double> r_mem;
  r_mem.set_n_mem(2 * (1 + par.size()), n_threads);

  // compute
  auto all_idx = get_indices(indices, *terms_ptr);
  int const * idx = &all_idx[0];
  int n_fails(0);

  openmp_exception_ptr exception_handler;
  pedmod::cdf_methods const meth = pedmod::get_cdf_methods(method);

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
{
#endif
  double * wmem    = r_mem.get_mem(),
         * var_est = wmem + 1 + par.size();
  std::fill(wmem   , wmem    + 1 + par.size(), 0.);
  std::fill(var_est, var_est + 1 + par.size(), 0.);

#ifdef _OPENMP
#pragma omp for schedule(static) reduction(+:n_fails)
#endif
  for(int i = 0; i < all_idx.size(); ++i)
    exception_handler.run([&]() -> void {
      if(idx[i] >= static_cast<int>(terms.size()))
        return;
      bool did_fail(false);
      double const w_i = has_weights ? c_weights[idx[i]] : 1;
      if(std::abs(w_i) < std::numeric_limits<double>::epsilon())
        return;

      *wmem += terms.at(idx[i]).gr(
        &par[0], wmem + 1, var_est, maxvls, abs_eps, rel_eps, minvls,
        do_reorder, use_aprx, did_fail, w_i, meth);
      n_fails += did_fail;
    });
#ifdef _OPENMP
}
#endif

  exception_handler.rethrow_if_error();

  // aggregate the result
  Rcpp::NumericVector grad(n_fix + n_scales),
                   std_est(n_fix + n_scales + 1);
  double ll(0.);
  for(unsigned i = 0; i < n_threads; ++i){
    double *wmem = r_mem.get_mem(i);
    ll += *wmem;
    for(unsigned j = 0; j < par.size(); ++j){
      grad   [j] += wmem[j + 1];
      std_est[j] += wmem[j + 1 + par.size()];
    }
    std_est[par.size()] += wmem[par.size() + 1 + par.size()];
  }

  for(unsigned j = 0; j < par.size() + 1; ++j)
    std_est[j] = std::sqrt(std_est[j]);
  for(unsigned j = 1 + n_fix; j < par.size() + 1; ++j)
    std_est[j] *= par[j - 1];

  // account for exp(...)
  for(int i = n_fix; i < n_fix + n_scales; ++i)
    grad[i] *= par[i];
  grad.attr("logLik")  = Rcpp::NumericVector::create(ll);
  grad.attr("n_fails") = Rcpp::IntegerVector::create(n_fails);
  grad.attr("std")     = std_est;

  return grad;
}
