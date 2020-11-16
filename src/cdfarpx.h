#ifndef CDFAPRX_H
#define CDFAPRX_H

#include "arma-wrap.h"
#include <array>
#include <limits>
#include <memory>
#include <cmath>
#include "pnorm.h"
#include "qnorm.h"
#include "config.h"
#include <algorithm>
#include "new-mvt.h"
#include "norm-cdf-approx.h"
#include <cmath>

namespace pedmod {
extern "C"
{
  /**
   * @param N Dimension of the integral.
   * @param lower N-vector with lower bounds.
   * @param upper N-vector with upper bounds.
   * @param delta N-vector with mean.
   * @param correl N(N - 1)/2-dimensional vector with  upper triangle of the
   * correlation matrix.
   * @param infin N-dimensional vector indicating whether the bounds are
   * finite.
   * @param pivot not sure. Set it to true.
   * @param y N-dimensional vector with workig memory.
   * @param ND N unless there is double infinite regions.
   * @param A potentially permutated version of lower.
   * @param B potentially permutated version of upper.
   * @param DL potentially permutated version of delta.
   * @param cov N(N + 1)/2-dimensional vector with potentially permutated
   * Cholesky decomposition of correl.
   * @param infi potentially permutated version of infin.
   * @param inform non-zero if something went wrong.
   * @param idx N-dimensional vector with indices of applied permutation.
   * @param doscale logical for whether to scale the cholesky decomposition
   * to have ones in the diagonal.
   *
   * cov is scaled such that the diagonal entries are one. This implies that
   * it is __not__ the Cholesky decomposition of the correlation matrix.
   * A, B, and DL are scaled accordingly.
   */
  void F77_NAME(mvsort)(
      int const* /* N */, double const* /* lower */,
      double const* /* upper */, double const* /* delta */,
      double const* /* correl */, int const* /* infin */,
      double const* /* y */, int const* /* pivot */,
      int* /* ND */, double* /* A */, double* /* B */, double* /* DL */,
      double* /* cov */, int* /* infi */, int* /* inform */,
      int* /* idx */, int const* /* doscale */);
}

/**
 * @param lower The lower bounds.
 * @param upper The upper bounds.
 * @return the infin argument for the mvtdst subroutine.
 */
arma::ivec get_infin(arma::vec const &lower, arma::vec const &upper);

struct cor_vec_res {
  arma::vec cor_vec, sds;
};

/**
 returns the minum number of samples as the original Fortran code.
 */
inline int default_minvls(int dim){
  dim = std::max(1, dim);
  constexpr int const def_vals[10] =
    { 16L * 31L - 1L, 16L * 47L - 1L, 16L * 73L - 1L, 16L * 113L - 1L, 16L * 173L - 1L, 16L * 263L - 1L, 16L * 397L - 1L, 16L * 593L - 1L, 16L * 907L - 1L, 16L * 1361L - 1L };
  return def_vals[
    std::min(static_cast<int>(dim - 1L), 9)];
}

/**
 * @return a struct with the correlation matrix and standard deviation. The
 * correlation matrix is stored as a upper diagonal matrix.
 */
cor_vec_res get_cor_vec(const arma::mat&);

inline double safe_qnorm_w(double const x) noexcept {
  constexpr double const eps =
    std::numeric_limits<double>::epsilon() *
    std::numeric_limits<double>::epsilon() *
    std::numeric_limits<double>::epsilon();
  if(x <= 0)
    return  qnorm_w(eps, 0, 1, 1L, 0L);
  else if(x >= 1.)
    return -qnorm_w(eps, 0, 1, 1L, 0L);

  return qnorm_w   (x  , 0, 1, 1L, 0L);
}

inline double safe_qnorm_aprx(double const x) noexcept {
  constexpr double const eps =
    std::numeric_limits<double>::epsilon() *
    std::numeric_limits<double>::epsilon() *
    std::numeric_limits<double>::epsilon();
  if(x <= 0)
    return  qnorm_aprx(eps);
  else if(x >= 1.)
    return -qnorm_aprx(eps);

  return qnorm_aprx   (x);
}

/**
 * copies the upper triangular matrix.
 *
 * @param X Matrix top copy.
 * @param x Pointer to copy to.
 */
inline void copy_upper_tri
  (arma::mat const &X, double * __restrict__ x) noexcept {
  size_t const p = X.n_cols;
  for(unsigned c = 0; c < p; c++)
    for(unsigned r = 0; r <= c; r++, x++)
      *x = X.at(r, c);
}

/**
 * copies the lower triangular matrix.
 *
 * @param X Matrix top copy.
 * @param x Pointer to copy to.
 */
inline void copy_lower_tri
  (arma::mat const &X, double * __restrict__ x) noexcept {
  size_t const p = X.n_cols;
  for(unsigned c = 0; c < p; c++)
    for(unsigned r = c; r < p; r++, x++)
      *x = X.at(r, c);
}

/**
 * TODO: describe what this class does.
 */
template<class T_Functor, class out_type = typename T_Functor::out_type>
class cdf {
  T_Functor &functor;
  const int ndim,
            n_integrands;
  const bool use_aprx;
  bool is_permutated = false;

  static constexpr bool const
    needs_last_unif = T_Functor::needs_last_unif();

  // TODO: memory allocations
  arma::ivec infin;
  arma::ivec indices;
  arma::vec lower = arma::vec(ndim),
            upper = arma::vec(ndim),
       sigma_chol = arma::vec((ndim * (ndim + 1L)) / 2L),
             draw = arma::vec(ndim);

public:
  cdf(T_Functor &functor, arma::vec const &lower_in,
      arma::vec const &upper_in, arma::vec const &mu_in,
      arma::mat const &sigma_in, bool const do_reorder,
      bool const use_aprx):
    functor(functor),
    ndim(mu_in.n_elem),
    n_integrands(functor.get_n_integrands(mu_in, sigma_in)),
    use_aprx(use_aprx),
    infin(get_infin(lower_in, upper_in)),
    indices(ndim) {
    /* checks */
#ifdef DO_CHECKS
    if(sigma_in.n_cols != static_cast<size_t>(ndim) or
         sigma_in.n_rows != static_cast<size_t>(ndim))
      throw std::invalid_argument("cdf::cdf: invalid 'sigma_in'");
    if(n_integrands <= 0L)
      throw std::invalid_argument("cdf::cdf: invalid 'n_integrands'");
    if(lower_in.n_elem !=  static_cast<size_t>(ndim))
      throw std::invalid_argument("cdf::cdf: invalid 'lower_in'");
    if(upper_in.n_elem !=  static_cast<size_t>(ndim))
      throw std::invalid_argument("cdf::cdf: invalid 'upper_in'");
#endif

    /* re-scale */
    arma::vec sds = arma::sqrt(arma::diagvec(sigma_in)),
               mu = mu_in / sds;

    lower  = lower_in;
    lower /= sds;
    lower -= mu;

    upper  = upper_in;
    upper /= sds;
    upper -= mu;

    is_permutated = false;
    {
      int *idx = indices.begin();
      for(int i = 0; i < ndim; ++i, ++idx)
        *idx = i;
    }

    if(do_reorder and ndim > 1L){
      std::unique_ptr<double[]> tmp_mem(new double[2 * ndim]);

      double * const y     = draw.begin(),
             * const A     = tmp_mem.get(),
             * const B     = A + ndim,
             * const DL    = sds.memptr(),
             * const delta = mu.begin();
      sds.zeros();

      auto const correl = get_cor_vec(sigma_in);
      int const pivot = 1L, doscale = 1L;
      int F_inform = 0L,
             nddim = ndim;
      std::fill(delta, delta + ndim, 0.);

      arma::ivec infi(ndim);

      F77_CALL(mvsort)(
        &ndim, lower.memptr(), upper.memptr(), delta,
        correl.cor_vec.memptr(), infin.begin(), y, &pivot, &nddim, A, B,
        DL, sigma_chol.memptr(), infi.memptr(), &F_inform, indices.begin(),
        &doscale);

      if(F_inform != 0)
        throw std::runtime_error("cdf::cdf: error in mvsort");

      {
        int const *prev = indices.begin();
        for(int i = 0; i < ndim; ++prev, ++i){
          if(*prev != i){
            is_permutated = true;
            break;
          }
        }
      }

      if(is_permutated){
        for(int i = 0; i < ndim; ++i){
          lower[i] = *(A + i);
          upper[i] = *(B + i);
          infin[i] = infi[i];
        }

        arma::uvec uidx(ndim);
        arma::vec mu_permu(ndim);
        for(int i = 0; i < ndim; ++i){
          uidx[i]     = indices[i];
          mu_permu[i] = mu_in[uidx[i]];
        }

        return;

      } else {
        for(int i = 0; i < ndim; ++i){
          lower[i] = *(A + i);
          upper[i] = *(B + i);
        }

      }

    } else if(!do_reorder and ndim > 1L) {
      arma::mat tmp(ndim, ndim);
      tmp = sigma_in;
      tmp.each_row() /= sds.t();
      tmp.each_col() /= sds;
      if(!arma::chol(tmp, tmp))
        sigma_chol.fill(std::numeric_limits<double>::infinity());
      else
        copy_upper_tri(tmp, sigma_chol.memptr());

      if(ndim > 1L){
        /* rescale such that choleksy decomposition has ones in the diagonal */
        double * sc = sigma_chol.begin();
        for(int i = 0; i < ndim; ++i){
          double const scal = *(sc + i);
          lower[i] /= scal;
          upper[i] /= scal;
          double * const sc_end = sc + i + 1L;
          for(; sc != sc_end; ++sc)
            *sc /= scal;
        }
      }
    }
  }

  /**
   TODO: add description.
   */
  void operator()(
      int const *ndim_in, double const * unifs, int const *n_integrands_in,
      double * __restrict__ integrand_val) PEDMOD_NOEXCEPT {
#ifdef DO_CHECKS
    if(*ndim_in         != ndim)
      throw std::invalid_argument("cdf::eval_integrand: invalid 'ndim_in'");
    if(*n_integrands_in != n_integrands)
      throw std::invalid_argument("cdf::eval_integrand: invalid 'n_integrands_in'");
#endif

    double * const __restrict__ out = integrand_val,
           * const __restrict__ dr  = draw.begin();

    double w(1.);
    double const * __restrict__ sc   = sigma_chol.begin(),
                 * __restrict__ lw   = lower.begin(),
                 * __restrict__ up   = upper.begin(),
                 * __restrict__ unif = unifs;
    int const *infin_j = infin.begin();
    /* loop over variables and transform them to truncated normal
     * variables */
    for(int j = 0; j < ndim; ++j, ++sc, ++lw, ++up, ++infin_j, ++unif){
      double su(0.);
      double const *d = dr;
      for(int i = 0; i < j; ++i, sc++, d++)
        su += *sc * *d;

      auto pnorm_use = [&](double const x){
        return use_aprx ? pnorm_approx(x) : pnorm_std(x, 1L, 0L);
      };
      double lim_l(0.),
             lim_u(1.);
      if(*infin_j == 0L)
        lim_u = pnorm_use(*up - su);
      else if(*infin_j == 1L)
        lim_l = pnorm_use(*lw - su);
      else if(*infin_j == 2L){
        lim_l = pnorm_use(*lw - su);
        lim_u = pnorm_use(*up - su);

      }

      if(lim_l < lim_u){
        w *= lim_u - lim_l;
        if(needs_last_unif or j + 1 < ndim){
          double const quant_val = lim_l + *unif * (lim_u - lim_l);
          *(dr + j) =
            use_aprx ?
            safe_qnorm_aprx(quant_val) :
            safe_qnorm_w   (quant_val);
        }

      } else {
        w = 0.;
        std::fill(dr + j, dr + ndim, 0.);
        break;

      }
    }

    /* evaluate the integrand and weigth the result. */
    functor(dr, out, indices.begin());

    double * o = out;
    for(int i = 0; i < n_integrands; ++i, ++o)
      *o *= w;
  }

  /**
   * Performs the approximation.
   *
   * @param maxvls Maximum number of function evaluations allowed.
   * @param abs_eps Required absolute accuracy.
   * @param rel_eps Required relative accuracy.
   * @param minvls Minimum number of samples.
   */
  out_type approximate
  (int const maxvls, double const abs_eps, double const rel_eps,
   int const minvls = 0L){
#ifdef DO_CHECKS
    if(abs_eps <= 0 and rel_eps <= 0)
      throw std::invalid_argument("cdf::approximate: no valid convergence threshold");
    if(maxvls <= 0L)
      throw std::invalid_argument("cdf::approximate: invalid 'maxvls'");
#endif

    // setup
    std::unique_ptr<double[]> int_apprx(new double[n_integrands]);
    auto sampler = parallelrng::get_unif_drawer();

    if(ndim == 1L){
      /* handle the one-dimensional case as a special case */
      functor.univariate(int_apprx.get(), lower[0], upper[0]);

      return functor.get_output(int_apprx.get(), 0, 0, 0);

    } else if(std::isinf(*sigma_chol.begin()))
      throw std::runtime_error("std::isinf(*sigma_chol.begin())");

    /* perform the approximation */
    auto res = rand_Korobov(
      *this, ndim, minvls, maxvls, n_integrands, abs_eps, rel_eps,
      int_apprx.get(), sampler);

    return functor.get_output(int_apprx.get(), res.minvls, res.inform,
                              res.abserr);
  }
};

/**
 * functor classes used as template argument for cdf used to approximate the
 * likelihood. */
class likelihood {
public:
  constexpr static int get_n_integrands
    (arma::vec const&, arma::mat const&) {
    return 1L;
  }

  inline void operator()
    (double const *, double * out, int const *)
    PEDMOD_NOEXCEPT {
#ifdef DO_CHECKS
    if(!out)
      throw std::invalid_argument("likelihood::operator(): invalid out");
#endif
    *out = 1;
  }

  constexpr static bool needs_last_unif() {
    return false;
  }

  inline static void univariate(double * out,
                                double const lw, double const ub) {
    double const p_ub = std::isinf(ub) ? 1 : pnorm_std(ub, 1L, 0L),
                 p_lb = std::isinf(lw) ? 0 : pnorm_std(lw, 1L, 0L);
    *out = p_ub - p_lb;
  }

  struct out_type {
    /**
     * minvls Actual number of function evaluations used.
     * inform INFORM = 0 for normal exit, when
     *             ABSERR <= MAX(ABSEPS, RELEPS*||finest||)
     *          and
     *             INTVLS <= MAXCLS.
     *        INFORM = 1 If MAXVLS was too small to obtain the required
     *        accuracy. In this case a value finest is returned with
     *        estimated absolute accuracy ABSERR. */
    int minvls, inform;
    /// Maximum norm of estimated absolute accuracy of finest
    double abserr;
    /// The likelihood approximation
    double likelihood;
  };

  out_type get_output(double const * res, int const minvls,
                      int const inform, double const abserr){
    return out_type { minvls, inform, abserr, *res };
  }
};
} // namespace pedmod

#endif
