#ifndef PEDIGREE_LL_H
#define PEDIGREE_LL_H

#include "cdfaprx.h"
#include <stdexcept>
#include <limits>

namespace pedmod {
class pedigree_ll_term {
  /// design matrix for the fixed effects
  arma::mat const X;

  static cache_mem<double> dmem;

  unsigned const max_n_sequences;

public:
  /// object to compute the likelihood factor
  pedigree_l_factor l_factor;

  unsigned const n_members = X.n_rows,
              n_fix_effect = X.n_cols;

  pedigree_ll_term(arma::mat const &X_in, arma::vec const &y,
                   std::vector<arma::mat> const &scale_mats,
                   unsigned const max_threads,
                   unsigned const max_n_sequences):
    X(([&](){
      if(X_in.n_rows != y.n_elem)
        throw std::invalid_argument("pedigree_ll_term::pedigree_ll_term: y and X's dimension do not match");

      arma::mat out = X_in;
      for(arma::uword i = 0; i < X_in.n_rows; ++i)
        if(y[i] > 0)
          out.row(i) *= -1;

      return out;
    })()),
    max_n_sequences(max_n_sequences),
    l_factor(([&](){
      // set cache
      cdf<likelihood       >::alloc_mem(y.n_elem, max_threads);
      cdf<pedigree_l_factor>::alloc_mem(y.n_elem, max_threads);

      arma::vec z(y.n_elem);
      for(arma::uword i = 0; i < y.n_elem; ++i)
        z[i] = y[i] > 0 ? 1 : -1;

      std::vector<arma::mat> out = scale_mats;
      for(auto &x : out){
        if(x.n_rows != z.n_elem or x.n_cols != z.n_elem)
          throw std::invalid_argument("pedigree_ll_term::pedigree_ll_term: invalid scale matrices");

        for(arma::uword j = 0; j < z.n_elem; ++j)
          for(arma::uword i = 0; i < j; ++i){
            double const mult = z[i] * z[j];
            x.at(i, j) *= mult;
            x.at(j, i) *= mult;
          }
      }

      return pedigree_l_factor(out, max_threads, X.t(), max_n_sequences);
    })()) {
    // checks
    if(l_factor.n_mem != n_members)
      throw std::invalid_argument("pedigree_ll_term::pedigree_ll_term: X and scale matrices dimension do not match");

    // set working memory
    likelihood::alloc_mem(n_members, max_threads, max_n_sequences);
    dmem.set_n_mem(
      3 * n_members + n_members * n_members, max_threads);
  }

  struct fn_res {
    double log_likelihood;
    double estimator_var;
  };

  /**
   * Approximates the log-likelihood term.
   */
  fn_res fn
    (double const * par, unsigned const maxvls, double const abs_eps,
     double const rel_eps, int minvls, bool const do_reorder,
     bool const use_aprx, bool &did_fail, cdf_methods const method){
    did_fail = true;
    arma::vec mu(dmem.get_mem(), n_members, false),
           lower(mu.end()      , n_members, false),
           upper(lower.end()   , n_members, false);
    std::fill(lower.begin(), lower.end(),
              -std::numeric_limits<double>::infinity());
    upper.zeros();

    arma::vec beta(const_cast<double *>(par), n_fix_effect, false);
    for(unsigned i = 0; i < n_members; ++i)
      mu[i] = arma::dot(beta, X.row(i));

    arma::mat sig(upper.end(), n_members, n_members, false);
    l_factor.setup(sig, par + n_fix_effect, 1., true);

    likelihood func;
    if(minvls < 0)
      minvls = std::min<int>(1000, 100 * n_members);
    auto const res = cdf<likelihood>(
      func, lower, upper, mu, sig, do_reorder, use_aprx).approximate(
          maxvls, abs_eps, rel_eps, method, minvls, max_n_sequences);

    did_fail = res.inform > 0;

    double const log_likelihood = std::log(res.likelihood);

    // crude variance estimator based on the delta rule
    double const sig_est = res.abserr * 2 / 7 / res.likelihood,
                 var_est = sig_est * sig_est;

    return { log_likelihood, var_est };
  }

  /**
   * Approximates the log-likelihood term and the derivative.
   */
  double gr
    (double const * par, double * d_par, double * var_est, unsigned const maxvls,
     double const abs_eps, double const rel_eps, int minvls,
     bool const do_reorder, bool const use_aprx, bool &did_fail,
     double const weight, cdf_methods const method){
    did_fail = true;
    arma::vec mu(dmem.get_mem(), n_members, false),
           lower(mu.end()      , n_members, false),
           upper(lower.end()   , n_members, false);
    std::fill(lower.begin(), lower.end(),
              -std::numeric_limits<double>::infinity());
    upper.zeros();

    arma::vec beta(const_cast<double *>(par), n_fix_effect, false);
    for(unsigned i = 0; i < n_members; ++i)
      mu[i] = arma::dot(beta, X.row(i));

    arma::mat sig(upper.end(), n_members, n_members, false);
    {
      l_factor.setup(sig, par + n_fix_effect, 1., true);
      pedmod::likelihood lfunc;
      auto const norm_const = pedmod::cdf<pedmod::likelihood>(
        lfunc, lower, upper, mu, sig, do_reorder, use_aprx).approximate(
            maxvls, abs_eps, std::min(1., 10. * rel_eps), method, minvls,
            max_n_sequences);

      l_factor.setup(sig, par + n_fix_effect, norm_const.likelihood, false);
    }

    if(minvls < 0)
      minvls = std::min<unsigned>(1000, 100 * n_members);
    auto const res = cdf<pedigree_l_factor>(
      l_factor, lower, upper, mu, sig, do_reorder, use_aprx).approximate(
          maxvls, abs_eps, rel_eps, method, minvls,
          max_n_sequences);

    // derivatives for the slopes and the scale parameters
    int const n_scales = l_factor.scale_mats.size();
    for(unsigned i = 0; i < n_fix_effect + n_scales; ++i)
      d_par[i] += weight * res.derivs[i];

    // add variance terms to var_est. The first one for the log likelihood is an
    // application of the delta rule
    var_est[0] += weight * weight * res.sd_errs[0] * res.sd_errs[0] /
      (res.likelihood * res.likelihood);
    for(unsigned i = 1; i < n_fix_effect + n_scales + 1; ++i)
      var_est[i] += weight * weight * res.sd_errs[i] * res.sd_errs[i];

    did_fail = res.inform > 0;
    return weight * std::log(res.likelihood);
  }
};

} // namespace pedmod

#endif
