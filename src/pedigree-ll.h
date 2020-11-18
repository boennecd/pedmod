#ifndef PEDIGREE_LL_H
#define PEDIGREE_LL_H

#include "cdfarpx.h"
#include <stdexcept>
#include <limits>

namespace pedmod {
class pedigree_ll_term {
  /// design matrix for the fixed effects
  arma::mat const X;

public:
  /// object to compute the likelihood factor
  pedigree_l_factor l_factor;

  int const n_members = X.n_rows,
         n_fix_effect = X.n_cols;

  pedigree_ll_term(arma::mat const &X_in, arma::vec const &y,
                   std::vector<arma::mat> const &scale_mats):
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

      return pedigree_l_factor(out);
    })()) {
    // checks
    if(l_factor.n_mem != n_members)
      throw std::invalid_argument("pedigree_ll_term::pedigree_ll_term: X and scale matrices dimension do not match");
  }

  /**
   * Approximates the log-likelihood term.
   */
  double fn
    (double const * par, int const maxvls, double const abs_eps,
     double const rel_eps, int minvls, bool const do_reorder,
     bool const use_aprx){
    // TODO: memory allocations
    arma::vec mu(n_members),
           lower(n_members),
           upper(n_members, arma::fill::zeros);
    std::fill(lower.begin(), lower.end(),
              -std::numeric_limits<double>::infinity());

    arma::vec beta(const_cast<double *>(par), n_fix_effect, false);
    for(int i = 0; i < n_members; ++i)
      mu[i] = arma::dot(beta, X.row(i));

    // TODO: memory allocations
    arma::mat sig(n_members, n_members);
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
    // TODO: memory allocations
    arma::vec mu(n_members),
    lower(n_members),
    upper(n_members, arma::fill::zeros);
    std::fill(lower.begin(), lower.end(),
              -std::numeric_limits<double>::infinity());

    arma::vec beta(const_cast<double *>(par), n_fix_effect, false);
    for(int i = 0; i < n_members; ++i)
      mu[i] = arma::dot(beta, X.row(i));

    // TODO: memory allocations
    arma::mat sig(n_members, n_members);
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
