#ifndef PEDIGREE_LL_H
#define PEDIGREE_LL_H

#include "cdfarpx.h"
#include <stdexcept>
#include <limits>

namespace pedmod {
class pedigree_ll_term {
  /// design matrix for the fixed effects
  arma::mat const X;

  static cache_mem<double> dmem;

public:
  /// object to compute the likelihood factor
  pedigree_l_factor l_factor;

  int const n_members = X.n_rows,
         n_fix_effect = X.n_cols;

  pedigree_ll_term(arma::mat const &X_in, arma::vec const &y,
                   std::vector<arma::mat> const &scale_mats,
                   unsigned const max_threads):
    X(([&](){
      if(X_in.n_rows != y.n_elem)
        throw std::invalid_argument("pedigree_ll_term::pedigree_ll_term: y and X's dimension do not match");

      arma::mat out = X_in;
      for(arma::uword i = 0; i < X_in.n_rows; ++i)
        if(y[i] > 0)
          out.row(i) *= -1;

      return out;
    })()),
    l_factor(([&](){
      // set cache
      cdf<likelihood       >::set_cache(y.n_elem, max_threads);
      cdf<pedigree_l_factor>::set_cache(y.n_elem, max_threads);

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

      return pedigree_l_factor(out, max_threads);
    })()) {
    // checks
    if(l_factor.n_mem != n_members)
      throw std::invalid_argument("pedigree_ll_term::pedigree_ll_term: X and scale matrices dimension do not match");

    // set working memory
    likelihood::set_cache(n_members, max_threads);
    dmem.set_n_mem(
      3 * n_members + n_members * n_members, max_threads);
  }

  /**
   * Approximates the log-likelihood term.
   */
  double fn
    (double const * par, int const maxvls, double const abs_eps,
     double const rel_eps, int minvls, bool const do_reorder,
     bool const use_aprx){
    arma::vec mu(dmem.get_mem(), n_members, false),
           lower(mu.end()      , n_members, false),
           upper(lower.end()   , n_members, false);
    std::fill(lower.begin(), lower.end(),
              -std::numeric_limits<double>::infinity());
    upper.zeros();

    arma::vec beta(const_cast<double *>(par), n_fix_effect, false);
    for(int i = 0; i < n_members; ++i)
      mu[i] = arma::dot(beta, X.row(i));

    arma::mat sig(upper.end(), n_members, n_members, false);
    l_factor.setup(sig, par + n_fix_effect, true);

    likelihood func;
    if(minvls < 0)
      minvls = std::min(1000, 100 * n_members);
    auto const res = cdf<likelihood>(
      func, lower, upper, mu, sig, do_reorder, use_aprx).approximate(
          maxvls, abs_eps, rel_eps, minvls);

    return std::log(res.likelihood);
  }

  /**
   * Approximates the log-likelihood term and the derivative.
   */
  double gr
    (double const * par, double * d_par, int const maxvls,
     double const abs_eps, double const rel_eps, int minvls,
     bool const do_reorder, bool const use_aprx){
    arma::vec mu(dmem.get_mem(), n_members, false),
           lower(mu.end()      , n_members, false),
           upper(lower.end()   , n_members, false);
    std::fill(lower.begin(), lower.end(),
              -std::numeric_limits<double>::infinity());
    upper.zeros();

    arma::vec beta(const_cast<double *>(par), n_fix_effect, false);
    for(int i = 0; i < n_members; ++i)
      mu[i] = arma::dot(beta, X.row(i));

    arma::mat sig(upper.end(), n_members, n_members, false);
    l_factor.setup(sig, par + n_fix_effect, false);

    if(minvls < 0)
      minvls = std::min(1000, 100 * n_members);
    auto const res = cdf<pedigree_l_factor>(
      l_factor, lower, upper, mu, sig, do_reorder, use_aprx).approximate(
          maxvls, abs_eps, rel_eps, minvls);

    // derivatives for the slopes
    for(int i = 0; i < n_fix_effect; ++i){
      double term(0.);
      for(int j = 0; j < n_members; ++j)
        term += res.derivs[j] * X.at(j, i);
      d_par[i] += term / res.likelihood;
    }

    // derivatives for the scale parameters
    double * rhs = d_par + n_fix_effect;
    double const * lhs = res.derivs.begin() + n_members;
    int const n_scales = l_factor.scale_mats.size();
    for(int i = 0; i < n_scales; ++i)
      *rhs++ += *lhs++ / res.likelihood;

    return std::log(res.likelihood);
  }
};

} // namespace pedmod

#endif
