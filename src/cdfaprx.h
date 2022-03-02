#ifndef CDFAPRX_H
#define CDFAPRX_H

#include "arma-wrap.h"
#include <array>
#include <vector>
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
#include <stdexcept>
#include "ped-mem.h"
#include "string"
#include "config.h"
#include "find-tilting-param.h"
#include <limits.h>
#include "qtnorm.h"

#include <R_ext/RS.h> // F77_NAME and F77_CALL

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

  void F77_NAME(dtpsv)
    (const char * /* uplo */, const char * /* trans */, const char * /* diag */,
     const int * /* n */, const double * /* ap */, double * /* x */,
     const int* /* incx */, size_t, size_t, size_t);

  void F77_NAME(dpptri)
    (const char * /* uplo */, const int * /* n */, double *ap,
     int * /* info */, size_t);
}

/**
 * @parma out Vector with result.
 * @param lower The lower bounds.
 * @param upper The upper bounds.
 */
arma::ivec get_infin
(arma::ivec &out, arma::vec const &lower, arma::vec const &upper);

struct cor_vec_res {
  arma::vec cor_vec, sds;
};

/**
 returns the minimum number of samples as the original Fortran code.
 */
inline int default_minvls(int dim){
  dim = std::max(1, dim);
  constexpr int def_vals[10] =
    { 16L * 31L - 1L, 16L * 47L - 1L, 16L * 73L - 1L, 16L * 113L - 1L, 16L * 173L - 1L, 16L * 263L - 1L, 16L * 397L - 1L, 16L * 593L - 1L, 16L * 907L - 1L, 16L * 1361L - 1L };
  return def_vals[
    std::min(static_cast<int>(dim - 1L), 9)];
}

/**
 * @return a struct with the correlation matrix and standard deviation. The
 * correlation matrix is stored as a upper diagonal matrix.
 */
cor_vec_res get_cor_vec(const arma::mat&);

template<bool lower_tail, bool use_log, bool use_aprx>
double pnorm_use
  (double const x) {
  if constexpr(use_aprx){
    double const pnrm{pnorm_approx(x)};
    if(use_log)
      return lower_tail ? log(pnrm) :  log1p(-pnrm);
    return lower_tail ? pnrm : 1 - pnrm;
  }
  return pnorm_std(x, lower_tail, use_log);
}

/**
 * copies the upper triangular matrix.
 *
 * @param X Matrix to copy.
 * @param x Pointer to copy to.
 */
inline void copy_upper_tri
  (arma::mat const &X, double * PEDMOD_RESTRICT x) noexcept {
  arma::uword const p = X.n_cols;
  for(arma::uword c = 0; c < p; c++)
    for(arma::uword r = 0; r <= c; r++, x++)
      *x = X.at(r, c);
}

/**
 * copies the lower triangular matrix.
 *
 * @param X Matrix to copy.
 * @param x Pointer to copy to.
 */
inline void copy_lower_tri
  (arma::mat const &X, double * PEDMOD_RESTRICT x) noexcept {
  arma::uword const p = X.n_cols;
  for(arma::uword c = 0; c < p; c++)
    for(arma::uword r = c; r < p; r++, x++)
      *x = X.at(r, c);
}

enum cdf_methods : int {
  Korobov = 0,
  Sobol = 1
};

inline cdf_methods get_cdf_methods(int const x){
  if(x < 0 or x > 1)
    throw std::invalid_argument("cdf_methods is not implemented");
  return static_cast<cdf_methods>(x);
}

/**
 * TODO: describe what this class does.
 */
template<class T_Functor, class out_type = typename T_Functor::out_type>
class cdf {
  T_Functor &functor;
  const arma::uword ndim,
                    n_integrands;
  const bool use_aprx;
  bool is_permutated = false;
  bool use_tilting;

  static constexpr bool
    needs_last_unif = T_Functor::needs_last_unif();

  // cached memory to use
  static cache_mem<int   > imem;
  static cache_mem<double> dmem;

  arma::ivec infin;
  arma::ivec indices;

  double * PEDMOD_RESTRICT const lower      = dmem.get_mem(),
         * PEDMOD_RESTRICT const upper      = lower + ndim,
         * PEDMOD_RESTRICT const sigma_chol = upper + ndim,
         * PEDMOD_RESTRICT const tilt_param =
         sigma_chol + (ndim * (ndim + 1L)) / 2L,
         * PEDMOD_RESTRICT const draw       = tilt_param + ndim,
         * PEDMOD_RESTRICT const dtmp_mem   = draw + ndim * n_qmc_seqs();

  // memory that can be used
  int * const itmp_mem{indices.end()};

  /**
   computes multiple integrands simultaneously.
   */
  template<bool with_tilting, bool with_aprx>
  void evaluate_intrands(
      unsigned const *ndim_in, double const * unifs,
      unsigned const *n_integrands_in,
      double * PEDMOD_RESTRICT integrand_val, unsigned const n_draws) PEDMOD_NOEXCEPT {
#ifdef DO_CHECKS
    if(*ndim_in         != ndim)
      throw std::invalid_argument("cdf::eval_integrand: invalid 'ndim_in'");
    if(*n_integrands_in != n_integrands)
      throw std::invalid_argument("cdf::eval_integrand: invalid 'n_integrands_in'");
#endif

    constexpr double Inf = std::numeric_limits<double>::infinity();

    double * const PEDMOD_RESTRICT out      = integrand_val,
           * const PEDMOD_RESTRICT dr       = draw,
           * const PEDMOD_RESTRICT su       = dtmp_mem,
           * const PEDMOD_RESTRICT w        = su + n_draws,
           * const PEDMOD_RESTRICT lim_l    = w + n_draws,
           * const PEDMOD_RESTRICT lim_u    = lim_l + n_draws,
           * const PEDMOD_RESTRICT lim_diff = lim_u + n_draws;

    if constexpr(with_tilting)
      std::fill(w, w + n_draws, 0);
    else
      std::fill(w, w + n_draws, 1);

    double const * PEDMOD_RESTRICT sc   = sigma_chol,
                 * PEDMOD_RESTRICT lw   = lower,
                 * PEDMOD_RESTRICT up   = upper;

    int const *infin_j = infin.begin();

    /* loop over variables and transform them to truncated normal
     * variables */
    for(arma::uword j = 0; j < ndim; ++j, ++sc, ++lw, ++up, ++infin_j){
      std::fill(su, su + n_draws, 0);
      {
        double * dri = dr;
        for(unsigned i = 0; i < j; ++i, ++sc)
          for(unsigned k = 0; k < n_draws; ++k, ++dri)
            su[k] += *sc * *dri;
      }

      if(*infin_j == 0L){
        std::fill(lim_l, lim_l + n_draws, with_tilting ? -Inf : 0);
        for(unsigned k = 0; k < n_draws; ++k)
          lim_u[k] = *up - su[k];

      } else if(*infin_j == 1L){
        std::fill(lim_u, lim_u  + n_draws,  with_tilting ? Inf : 1);
        for(unsigned k = 0; k < n_draws; ++k)
          lim_l[k] = *lw - su[k];

      } else
        for(unsigned k = 0; k < n_draws; ++k) {
          lim_l[k] = *lw - su[k];
          lim_u[k] = *up - su[k];
        }

      bool const is_last_run{j + 1 >= ndim};
      if constexpr(with_tilting){
        auto subtract_tilt = [&]{
          for(unsigned k = 0; k < n_draws; ++k) {
            lim_l[k] -= tilt_param[j];
            lim_u[k] -= tilt_param[j];
          }
        };

        if constexpr(needs_last_unif)
          subtract_tilt();
        else if(!is_last_run)
          subtract_tilt();
      } else {
        if(*infin_j == 0){
          for(unsigned k = 0; k < n_draws; ++k)
            lim_u[k] = pnorm_use<true, false, with_aprx>(lim_u[k]);

        } else if(*infin_j == 1){
          for(unsigned k = 0; k < n_draws; ++k)
            lim_l[k] = pnorm_use<true, false, with_aprx>(lim_l[k]);

        } else {
          for(unsigned k = 0; k < n_draws; ++k){
            lim_l[k] = pnorm_use<true, false, with_aprx>(lim_l[k]);
            lim_u[k] = pnorm_use<true, false, with_aprx>(lim_u[k]);
          }

        }
      }

      auto set_dr_n_weight = [&]{
        unsigned const offset = j * n_draws;

        if constexpr(with_tilting){
          for(unsigned k = 0; k < n_draws; ++k){
            double val, log_pnrms_diff;

            if(lim_l[k] > 0){
              double const v_lb{pnorm_use<false, true, with_aprx>(lim_l[k])},
                           v_ub{pnorm_use<false, true, with_aprx>(lim_u[k])};

              log_pnrms_diff = v_lb + std::log1p(-exp(v_ub - v_lb));

              val = std::exp(v_lb) -
                unifs[k * ndim + j] * std::exp(log_pnrms_diff);
              val = qnorm_w(val, 0, 1, 0, 0);

            } else if(lim_u[k] < 0){
              double const v_lb{pnorm_use<true, true, with_aprx>(lim_l[k])},
                           v_ub{pnorm_use<true, true, with_aprx>(lim_u[k])};

              log_pnrms_diff = v_ub + std::log1p(-exp(v_lb - v_ub));

              if(lim_u[k] < -35)
                val = qtnorm(unifs[k * ndim + j], lim_l[k], lim_u[k]);
              else {
                val =
                  std::exp(v_lb)
                    + unifs[k * ndim + j] * std::exp(log_pnrms_diff);
                val = -qnorm_w(val, 0, 1, 0, 0);
              }

            } else {
              double const v_lb{pnorm_use<true, false, with_aprx>(lim_l[k])},
                           v_ub{pnorm_use<false, false, with_aprx>(lim_u[k])};

              log_pnrms_diff = std::log1p(-v_lb - v_ub);

              val = pnorm_use<false, false, with_aprx>(lim_l[k]) -
                unifs[k * ndim + j] * std::exp(log_pnrms_diff);
              val = qnorm_w(val, 0, 1, 0, 0);

            }

            double const val_shifted{val + tilt_param[j]};
            dr[offset + k] = val_shifted;
            w[k] += log_pnrms_diff +
              (tilt_param[j] - 2 * val_shifted) * tilt_param[j] / 2;

          }
        } else {
          for(unsigned k = 0; k < n_draws; ++k)
            lim_diff[k] = lim_u[k] - lim_l[k];
          for(unsigned k = 0; k < n_draws; ++k)
            w[k] *= lim_diff[k];
          for(unsigned k = 0; k < n_draws; ++k){
            double const quant_val
              {lim_l[k] + unifs[k * ndim + j] * lim_diff[k]};
            dr[offset + k] = with_aprx ? qnorm_aprx(quant_val)
                                       : qnorm_w   (quant_val, 0, 1, 1, 0);
          }
        }

        for(unsigned k = 0; k < n_draws; ++k)
          if(lim_l[k] >= lim_u[k] ||
             unifs[k * ndim + j] <= 0 || unifs[k * ndim + j] >= 1){
            // for reasons that are beyond me at the moment, unifs are
            // (although very rarely) sometimes exactly 0 or 1 with some
            // methods which gives +/-Inf values when the appropriate limit
            // is +/-Inf
            w[k] = with_tilting ? -Inf : 0;
            dr[offset + k] = 0;
          }
      };

      if constexpr(needs_last_unif)
        set_dr_n_weight();
      else if(!is_last_run)
        set_dr_n_weight();
      else {
        for(unsigned k = 0; k < n_draws; ++k){
          if constexpr(with_tilting){
            double log_pnrms_diff;

            if(lim_l[k] > 0){
              double const v_lb{pnorm_use<false, true, with_aprx>(lim_l[k])},
                           v_ub{pnorm_use<false, true, with_aprx>(lim_u[k])};

              log_pnrms_diff = v_lb + std::log1p(-exp(v_ub - v_lb));

            } else if(lim_u[k] < 0){
              double const v_lb{pnorm_use<true, true, with_aprx>(lim_l[k])},
                           v_ub{pnorm_use<true, true, with_aprx>(lim_u[k])};

              log_pnrms_diff = v_ub + std::log1p(-exp(v_lb - v_ub));

            } else {
              double const v_lb{pnorm_use<true, false, with_aprx>(lim_l[k])},
                           v_ub{pnorm_use<false, false, with_aprx>(lim_u[k])};

              log_pnrms_diff = std::log1p(-v_lb - v_ub);

            }

            w[k] += log_pnrms_diff;

          } else
            w[k] *= lim_u[k] - lim_l[k];

        }

        for(unsigned k = 0; k < n_draws; ++k)
          if(lim_l[k] >= lim_u[k])
            w[k] = with_tilting ? -Inf : 0;
      }
    }

    /* evaluate the integrand and weight the result. */
    functor(dr, out, indices.begin(), is_permutated, n_draws);

    // multiply by the density
    double *o{out};
    for(unsigned k = 0; k < n_draws; ++k){
      if constexpr(with_tilting)
        w[k] = std::exp(w[k]);
      else if(std::isnan(w[k]))
        // the method is not numerically stable in very extreme settings
        // and we may have set a draw equal to +/-Inf which will give
        // a NaN
        w[k] = 0;

      w[k] /= functor.get_norm_constant();
      if(w[k] == 0){
        std::fill(o, o + n_integrands, 0);
        o += n_integrands;

      } else
        for(unsigned i = 0; i < n_integrands; ++i, ++o)
          *o *= w[k];
    }
  }

public:
  /**
   * must be called prior to calling the constructor or any member
   * functions.
   */
  static void alloc_mem(unsigned const max_ndim, unsigned const max_threads) {
    unsigned const n_up_tri = (max_ndim * (max_ndim + 1)) / 2;

    imem.set_n_mem(3 * max_ndim, max_threads);
    // TODO: this is wasteful. Look through this
    dmem.set_n_mem(
      (7 + n_qmc_seqs()) * max_ndim + n_up_tri + max_ndim * max_ndim +
        5 * n_qmc_seqs(),
      max_threads);
  }

  cdf(T_Functor &functor, arma::vec const &lower_in,
      arma::vec const &upper_in, arma::vec const &mu_in,
      arma::mat const &sigma_in, bool const do_reorder,
      bool const use_aprx, bool const use_tilting_in):
    functor(functor),
    ndim(mu_in.n_elem),
    n_integrands(functor.get_n_integrands()),
    use_aprx(use_aprx),
    use_tilting{use_tilting_in},
    infin(([&](){
      arma::ivec out(imem.get_mem(), ndim, false);
      get_infin(out, lower_in, upper_in);
      return out;
    })()),
    indices(infin.end(), ndim, false) {
    /* checks */
    if(lower_in.n_elem > 1000 or lower_in.n_elem < 1)
      throw std::invalid_argument("cdf<T_Functor, out_type>: Either dimension zero or dimension greater than 1000");

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
    double * cur_dtmp_mem = dtmp_mem;
    auto get_dmem = [&](unsigned const n_ele) -> double * {
      double * out = cur_dtmp_mem;
      cur_dtmp_mem += n_ele;
      return out;
    };

    double * sds = get_dmem(ndim);
    for(arma::uword i = 0; i < ndim; ++i){
      sds[i] = std::sqrt(sigma_in.at(i, i));
      lower[i] = (lower_in[i] - mu_in[i]) / sds[i];
      upper[i] = (upper_in[i] - mu_in[i]) / sds[i];
    }

    is_permutated = false;
    {
      int *idx = indices.begin();
      for(arma::uword i = 0; i < ndim; ++i, ++idx)
        *idx = static_cast<int>(i);
    }

    auto setup_titling = [&]{
      if(use_tilting){
        for(arma::uword i = 0; i < ndim; ++i){
          if(infin[i] == 0)
            lower[i] = -std::numeric_limits<double>::infinity();
          else if(infin[i] == 1)
            upper[i] = std::numeric_limits<double>::infinity();
        }

        if(ndim < 2)
          use_tilting = false;
        else {
          auto find_tilt_res = find_tilting_param
            (ndim, lower, upper, sigma_chol, 1e-8);

          use_tilting = find_tilt_res.success;
          if(find_tilt_res.success)
            std::copy
            (find_tilt_res.tilting_param.begin(),
             find_tilt_res.tilting_param.end(), tilt_param);
        }
      }
    };

    if(do_reorder and ndim > 1L){
      double * const y     = draw,
             * const A     = get_dmem(ndim),
             * const B     = get_dmem(ndim),
             * const DL    = sds,
             * const delta = get_dmem(ndim);
      std::fill(sds, sds + ndim, 0.);

      auto const correl = get_cor_vec(sigma_in);
      int const pivot = 1L, doscale = 1L;
      int F_inform = 0L,
             nddim = static_cast<int>(ndim);
      std::fill(delta, delta + ndim, 0.);
      arma::ivec infi(itmp_mem, ndim, false);

      {
        int ndim_int(static_cast<int>(ndim));
        F77_CALL(mvsort)(
          &ndim_int, lower, upper, delta,
          correl.cor_vec.memptr(), infin.begin(), y, &pivot, &nddim, A, B,
          DL, sigma_chol, infi.memptr(), &F_inform, indices.begin(),
          &doscale);
      }

      if(F_inform != 0)
        throw std::runtime_error("cdf::cdf: error in mvsort");

      for(arma::uword i = 0; i < ndim; ++i){
        if(indices[i] != static_cast<int>(i)){
          is_permutated = true;
          break;
        }
      }

      if(is_permutated){
        for(arma::uword i = 0; i < ndim; ++i){
          lower[i] = A[i];
          upper[i] = B[i];
          infin[i] = infi[i];
        }

        arma::mat sigma_permu(get_dmem(ndim * ndim), ndim, ndim, false);
        for(arma::uword j = 0; j < ndim; ++j)
          for(arma::uword i = 0; i < ndim; ++i)
            sigma_permu.at(i, j) = sigma_in.at(indices[i], indices[j]);

        setup_titling();
        functor.prep_permutated(sigma_permu, indices.begin());

        return;

      } else
        for(arma::uword i = 0; i < ndim; ++i){
          lower[i] = A[i];
          upper[i] = B[i];
        }

    } else if(!do_reorder and ndim > 1L) {
      arma::mat tmp(get_dmem(ndim * ndim), ndim, ndim, false);
      tmp = sigma_in;
      for(arma::uword i = 0; i < ndim; ++i)
        for(arma::uword j = 0; j < ndim; ++j)
          tmp.at(i, j) /= sds[i] * sds[j];

      if(!arma::chol(tmp, tmp)) // TODO: memory allocation
        std::fill(sigma_chol, sigma_chol + (ndim * (ndim + 1L)) / 2L,
                  std::numeric_limits<double>::infinity());
      else
        copy_upper_tri(tmp, sigma_chol);

      if(ndim > 1L){
        /* rescale such that Choleksy decomposition has ones in the diagonal */
        double * sc = sigma_chol;
        for(arma::uword i = 0; i < ndim; ++i){
          double const scal = sc[i];
          lower[i] /= scal;
          upper[i] /= scal;
          double * const sc_end = sc + i + 1L;
          for(; sc != sc_end; ++sc)
            *sc /= scal;
        }
      }
    } else
      *sigma_chol = 1.;

    setup_titling();
    functor.prep_permutated(sigma_in, indices.begin());
  }

  /**
   TODO: add description.
   */
  void operator()(
      unsigned const *ndim_in, double const * unifs,
      unsigned const *n_integrands_in,
      double * PEDMOD_RESTRICT integrand_val, unsigned const n_draws) PEDMOD_NOEXCEPT {
    if(use_tilting) {
      if(use_aprx)
        evaluate_intrands<true, true>
         (ndim_in, unifs, n_integrands_in, integrand_val,n_draws);
      else
        evaluate_intrands<true, false>
          (ndim_in, unifs, n_integrands_in, integrand_val,n_draws);
    } else {
      if(use_aprx)
        evaluate_intrands<false, true>
          (ndim_in, unifs, n_integrands_in, integrand_val,n_draws);
      else
        evaluate_intrands<false, false>
          (ndim_in, unifs, n_integrands_in, integrand_val,n_draws);
    }

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
  (size_t const maxvls, double const abs_eps, double const rel_eps,
   cdf_methods const method, size_t const minvls, unsigned const n_sequences){
#ifdef DO_CHECKS
    if(abs_eps <= 0 and rel_eps <= 0)
      throw std::invalid_argument("cdf::approximate: no valid convergence threshold");
    if(maxvls <= 0L)
      throw std::invalid_argument("cdf::approximate: invalid 'maxvls'");
#endif

    // setup
    // needs to have at least n_integrands memory to use.
    double * const int_apprx = functor.get_wk_mem(),
           * const int_sdest = int_apprx + n_integrands;

    auto sampler = parallelrng::get_unif_drawer();

    if(ndim == 1L){
      /* handle the one-dimensional case as a special case */
      functor.univariate(int_apprx, lower[0], upper[0]);
      indices[0] = 0;
      // assume that there is zero error in the univariate case
      std::fill(int_sdest, int_sdest + n_integrands, 0.);

      return functor.get_output(int_apprx, int_sdest, 0, 0, 0,
                                indices.begin());

    } else if(std::isinf(*sigma_chol))
      throw std::runtime_error("std::isinf(*sigma_chol.begin())");

    /* perform the approximation */
    auto res = ([&]() -> rand_Korobov_output {
      if(method == cdf_methods::Sobol)
        return sobol_wrapper<cdf<T_Functor> >::comp(
            *this, ndim, minvls, maxvls, n_integrands, abs_eps, rel_eps,
            int_apprx, int_sdest, sampler, sobol::scrambling_type::owen,
            n_sequences);
      if(method != cdf_methods::Korobov)
        throw std::invalid_argument("method is not implemented");

      return rand_Korobov<cdf<T_Functor> >::comp(
          *this, ndim, minvls, maxvls, n_integrands, abs_eps, rel_eps,
          int_apprx, int_sdest, sampler, n_sequences);
    })();

    return functor.get_output(int_apprx, int_sdest, res.minvls, res.inform,
                              res.abserr, indices.begin());
  }
};

template<class T_Functor, class out_type>
cache_mem<int   > cdf<T_Functor, out_type>::imem;
template<class T_Functor, class out_type>
cache_mem<double> cdf<T_Functor, out_type>::dmem;

/**
 * functor classes used as template argument for cdf used to approximate the
 * likelihood. */
class likelihood {
  static cache_mem<double> dmen;

public:
  static void alloc_mem
  (unsigned const max_dim, unsigned const max_threads,
   unsigned const max_n_sequences){
    rand_Korobov<cdf<likelihood> >::alloc_mem(
        max_dim, get_n_integrands(), max_threads);
    sobol_wrapper<cdf<likelihood> >::alloc_mem(
        max_dim, get_n_integrands(), max_threads, max_n_sequences);
    dmen.set_n_mem(2, max_threads);
  }

  double * get_wk_mem(){
    return dmen.get_mem();
  }

  constexpr static unsigned get_n_integrands() {
    return 1;
  }

  void operator()
    (double const *, double * out, int const *, bool const,
     unsigned const n_draws)
    PEDMOD_NOEXCEPT {
#ifdef DO_CHECKS
    if(!out)
      throw std::invalid_argument("likelihood::operator(): invalid out");
#endif
    std::fill(out, out + n_draws, 1);
  }

  constexpr static bool needs_last_unif() {
    return false;
  }

  constexpr static double get_norm_constant() {
    return 1;
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
    size_t minvls;
    int inform;
    /// maximum norm of estimated absolute accuracy of finest
    double abserr;
    /// likelihood approximation
    double likelihood;
  };

  out_type get_output(double const * res, double const * sdest,
                      size_t const minvls, int const inform, double const abserr,
                      int const *){
    return out_type { minvls, inform, abserr, *res };
  }

  void prep_permutated(arma::mat const&, int const*) { }
};

/**
 * functor classes used as template argument for cdf used to approximate the
 * derivatives of the likelihood factors for each family. That is, the
 * likelihood
 *
 *   int_B phi(u;X.beta, I + sum_k sigma[k] * C[k]) du
 *
 * For given matrices C and a design matrix X.
 *
 * The returned approximations is a) the likelihood factor, b) the
 * derivative of log likelihood w.r.t. the fixed effect coefficients and w.r.t.
 * each of the scale parameters.
 */
class pedigree_l_factor {
public:
  /// the scale matrices for the different effects
  std::vector<arma::mat> const scale_mats;
  /// the number of members in this family
  arma::uword const n_mem = scale_mats[0].n_rows;
  /// design matrix for the fixed effects
  arma::mat const X;
  /// the number of fixed effects
  arma::uword const n_fix = X.n_cols,
                 n_scales = scale_mats.size(),
             n_integrands = 1 + n_fix + n_scales;
  /// scale free constant to check that a matrix is positive semi definite
  static constexpr double eps_pos_def =
    10 * std::numeric_limits<double>::epsilon();

private:
  /// working memory
  static cache_mem<double> dmem;
  static cache_mem<int   > imem;

  /**
   * the matrix  [cluster size] x [number of fixed effect] matrix which is
   * needed to compute the derivatives w.r.t. the slopes for the fixed effects.
   */
  double * d_fix_mat;

  /**
   * Let S^top S = Sigma be the Cholesky decomposition and C_{i1}, ..., C_{iK}
   * be the scale matrices. Then this is the first pointer to the following for
   * each of the K scale matrices:
   *  a. the upper triangular part of the Cholesky decomposition of
   *     S^-top C_{i1}S^-1, ..., S^-top C_{iK}S^-1 if the matrix
   *     S^-top C_{ik}S^-1 is positive definite.
   *  b. the transpose of the Eigen vectors Q scaled by the square root of the
   *     Eigen values where QDQ^top = S^-top C_{i1}S^-1.
   *
   * These objects can be found increments of n_mem * n_mem
   */
  double * S_C_S_matrices;

  /**
   * stores how many non-zero Eigen values each S_C_S_matrices has. It is minus
   * one if it is a Cholesky decomposition is used.
   */
  int * S_C_n_eigen;

  /// points to the upper triangular part of the inverse.
  double * sig_inv;

  /// working memory to be used by cdf
  double * cdf_mem;

  /// working memory that can be used for anything
  double * interal_mem;

  /// array of pointer to the scale matrices' element which we will need.
  std::unique_ptr<double const *[]> scale_mats_ptr =
    std::unique_ptr<double const *[]>(new double const *[n_scales]);

  /// the normalization constant
  double norm_const = std::numeric_limits<double>::quiet_NaN();

public:
  /// sets the scale matrices. There are no checks on the validity
  pedigree_l_factor(std::vector<arma::mat> const &scale_mats,
                    unsigned const max_threads, arma::mat const &X_in,
                    unsigned const max_n_sequences):
  scale_mats(scale_mats), X(X_in.t()) {
    // checks
    if(scale_mats.size() < 1)
      throw std::invalid_argument("pedigree_l_factor::pedigree_l_factor: not scale matrices are passed");
    arma::uword const u_mem = n_mem;
    for(auto &S : scale_mats)
      if(S.n_rows != u_mem or S.n_rows != u_mem)
        throw std::invalid_argument("pedigree_l_factor::pedigree_l_factor: not all scale matrices are square matrices or have the same dimensions");
    if(X.n_rows != u_mem)
      throw std::invalid_argument("pedigree_l_factor::pedigree_l_factor: invalid X");

    // setup working memory
    rand_Korobov<cdf<pedigree_l_factor> >::alloc_mem(
        n_mem, get_n_integrands(), max_threads);
    sobol_wrapper<cdf<pedigree_l_factor> >::alloc_mem(
        n_mem, get_n_integrands(), max_threads, max_n_sequences);
    size_t const working_memory =
      std::max<size_t>(
        // for prep_permutated
        3 * n_mem * n_mem + n_mem + n_fix * n_mem,
        // for operator()
        2 * n_qmc_seqs());
    dmem.set_n_mem(
       (n_mem * (n_mem + 1)) / 2 +
        n_mem * n_mem * scale_mats.size() +
        n_fix * n_mem +
        2 * get_n_integrands() +
        working_memory,
      max_threads);
    imem.set_n_mem(n_scales, max_threads);

    // check that the scale matrices are positive semi definite
    arma::vec vdum(dmem.get_mem(), n_mem, false);
    for(arma::mat const &m : scale_mats){
      if(!arma::eig_sym(vdum, m))
        throw std::runtime_error("Eigen decomposition failed in pedigree_l_factor constructor");

      std::reverse(vdum.begin(), vdum.end());
      if(vdum[0] < 0)
        throw std::invalid_argument  ("None positive definite scale matrix. The largest eigen value is " +
                                      std::to_string(vdum[0]));
      double const eps = eps_pos_def * n_mem * vdum[0];
      for(arma::uword i = 1; i < n_mem; ++i)
        if(vdum[i] < -eps)
          throw std::invalid_argument(
              "None positive definite scale matrix. Largest eigen value is " +
                std::to_string(vdum[0]) + " and one eigen value is equal to " +
                std::to_string(vdum[i]));
    }
  }

  unsigned get_n_integrands() PEDMOD_NOEXCEPT {
    return n_integrands;
  }

  double * get_wk_mem() PEDMOD_NOEXCEPT {
    return cdf_mem;
  }

  constexpr static bool needs_last_unif() PEDMOD_NOEXCEPT {
    return true;
  }

  double get_norm_constant() PEDMOD_NOEXCEPT {
    return norm_const;
  }

  /**
   * setups the covariance matrix to use. This method must be called be
   * prior to using the object in an approximation.
   *
   * Args:
   *   sig: the covariance matrix.
   *   scales: scale parameters.
   *   norm_constant_arg: the normalization constant.
   *   only_cov: true if only the covariance matrix should be computed
   */
  void setup(arma::mat &sig, double const *scales,
             double const norm_constant_arg,
             bool const only_cov = false){
    sig.zeros(n_mem, n_mem);
    sig.diag() += 1;
    for(unsigned i = 0; i < scale_mats.size(); ++i)
      sig += scales[i] * scale_mats[i];
    if(only_cov)
      return;
    norm_const = norm_constant_arg;

    double * next_mem = dmem.get_mem();
    sig_inv = next_mem;
    next_mem += (n_mem * (n_mem + 1)) / 2;
    cdf_mem = next_mem;
    next_mem += 2 * get_n_integrands();

    arma::mat t1(next_mem, n_mem, n_mem, false);
    if(!arma::inv_sympd(t1, sig))
      throw std::runtime_error("pedigree_ll_factor::setup: inv_sympd failed");

    copy_upper_tri(t1, sig_inv);
  }

  void prep_permutated(arma::mat const &sig, int const *indices) {
    if(n_mem < 2)
      return;

    // create the objects we need. We have to account for the memory used by
    // setup
    double * next_mem = cdf_mem + 2 * get_n_integrands();
    interal_mem = next_mem + n_fix * n_mem + n_mem * n_mem * n_scales;

    arma::vec dum_vec(interal_mem  , n_mem, false);
    arma::mat t1(dum_vec.end(), n_mem, n_mem, false),
              t2(t1.end()     , n_mem, n_mem, false),
              t3(t2.end()     , n_mem, n_mem, false),
         X_permu(t3.end()     , n_mem, n_fix, false);
    if(!arma::chol(t1, sig, "lower"))
      throw std::runtime_error("pedigree_ll_factor::setup: chol failed");

    // permute X
    for(arma::uword j = 0; j < n_fix; ++j)
      for(arma::uword i = 0; i < n_mem; ++i)
        X_permu(i, j) = X(indices[i], j);

    d_fix_mat = next_mem;
    arma::mat d_fix_obj(d_fix_mat, n_mem, n_fix, false);
    next_mem += n_mem * n_fix;
    arma::solve(d_fix_obj, arma::trimatl(t1), X_permu);

    // set up the array we need
    S_C_n_eigen = imem.get_mem();
    S_C_S_matrices = next_mem;
    {
      unsigned s(0);
      for(arma::mat const &m : scale_mats){
        // setup the permutation of the scale matrix
        for(arma::uword j = 0; j < n_mem; ++j)
          for(arma::uword i = 0; i < n_mem; ++i)
            t2(i, j) = m(indices[i], indices[j]);

        arma::solve(t3, arma::trimatl(t1), t2);
        arma::inplace_trans(t3);
        arma::solve(t2, arma::trimatl(t1), t3);

        if(arma::chol(t3, t2)){
          S_C_n_eigen[s] = -1;
          arma::inplace_trans(t3);
          copy_lower_tri(t3, next_mem);

        } else {
          arma::mat &eigen_vectors = t3;

          // form the Eigen decomposition
          if(!arma::eig_sym(dum_vec, eigen_vectors, t2))
            throw std::runtime_error("Eigen decomposition failed");

          // count the number of Eigen values greater than zero
          double const eps = eps_pos_def * n_mem * dum_vec[n_mem - 1];
          arma::uword j = n_mem - 1;
          for(; j > 0; --j)
            if(dum_vec[j - 1] < eps)
              break;
          arma::uword const n_eigen_vectors = n_mem - j;
          S_C_n_eigen[s] = static_cast<int>(n_eigen_vectors);

          // over write the first column
          arma::uword const n_unused = n_mem - n_eigen_vectors;
          for(arma::uword k = 0; k < n_eigen_vectors; ++k){
            dum_vec[k] = dum_vec[k + n_unused];
            std::copy(eigen_vectors.colptr(k + n_unused),
                      eigen_vectors.colptr(k + n_unused) + n_mem,
                      eigen_vectors.colptr(k));
          }

          // scale the Eigen vectors
          double * eg_val = eigen_vectors.begin();
          for(arma::uword k = 0; k < n_eigen_vectors; ++k){
            double const scale = std::sqrt(dum_vec[k]);
            for(arma::uword j = 0; j < n_mem; ++j, ++eg_val)
              *eg_val *= scale;
          }

          std::copy(eigen_vectors.begin(),
                    eigen_vectors.begin() + n_mem * n_eigen_vectors, next_mem);
        }

        ++s;
        next_mem += n_mem * n_mem;
      }
    }
  }

  void operator()
    (double const * PEDMOD_RESTRICT draw, double * PEDMOD_RESTRICT out,
     int const *, bool const, unsigned const n_draws) {
      for(unsigned k = 0; k < n_draws; ++k)
        out[k * n_integrands] = 1;

      double * PEDMOD_RESTRICT sum       = interal_mem,
             * PEDMOD_RESTRICT sqrt_term = sum + n_draws;

      // derivatives w.r.t. the fixed effects
      {
        double const *xij = d_fix_mat;
        for(unsigned j = 0; j < n_fix; ++j){
          double const *d = draw;
          std::fill(sum, sum + n_draws, 0);

          for(arma::uword i = 0; i < n_mem; ++i, ++xij)
            for(unsigned k = 0; k < n_draws; ++k, ++d)
              sum[k] += *xij * *d;

          unsigned const offset = j + 1;
          for(unsigned k = 0; k < n_draws; ++k)
            out[offset + k * n_integrands] = sum[k];
        }
      }

      // handle the derivatives w.r.t. the scale parameters
      double * next_mat = S_C_S_matrices;
      for(arma::uword s = 0; s < n_scales; ++s, next_mat += n_mem * n_mem){
        std::fill(sum, sum + n_draws, 0);

        if(S_C_n_eigen[s] < 0){
          // a Cholesky decomposition was used
          double *chol_ele = next_mat;

          for(arma::uword c = 0; c < n_mem; ++c){
            std::fill(sqrt_term, sqrt_term + n_draws, 0);
            double const * d = draw + c * n_draws;
            for(arma::uword r = 0; r < n_mem - c; ++r, ++chol_ele)
              for(unsigned k = 0; k < n_draws; ++k, ++d)
                sqrt_term[k] += *chol_ele * *d;

            for(unsigned k = 0; k < n_draws; ++k)
              sum[k] += sqrt_term[k] * sqrt_term[k];
          }

        } else {
          // an Eigen decomposition was used
          double *eig_vec_ele = next_mat;
          unsigned const n_eigen_vec = static_cast<unsigned>(S_C_n_eigen[s]);

          for(unsigned r = 0; r < n_eigen_vec; ++r){
            std::fill(sqrt_term, sqrt_term + n_draws, 0);
            double const * d = draw;
            for(arma::uword c = 0; c < n_mem; ++c, ++eig_vec_ele)
              for(unsigned k = 0; k < n_draws; ++k, ++d)
                sqrt_term[k] += *eig_vec_ele * *d;

            for(unsigned k = 0; k < n_draws; ++k)
              sum[k] += sqrt_term[k] * sqrt_term[k];
          }
        }

        unsigned const offset = 1 + n_fix + s;
        for(unsigned k = 0; k < n_draws; ++k)
          out[offset + k * n_integrands] = sum[k] * .5;
      }
    }

  void univariate(double * out, double const lw, double const ub) {
    constexpr double log_sqrt_2_pi_inv{0.918938533204673};
    auto log_dnrm = [&](double const x){
      return -x * x / 2. - log_sqrt_2_pi_inv;
    };

    bool const f_ub = std::isinf(ub),
               f_lb = std::isinf(lw);

    double const p_ub = f_ub ? 1 : pnorm_std(ub, 1L, 0L),
                 p_lb = f_lb ? 0 : pnorm_std(lw, 1L, 0L),
                 d_ub = f_ub ? 0 : std::exp(log_dnrm(ub) - pnorm_std(ub, 1L, 1L)),
                 d_lb = f_lb ? 0 : std::exp(log_dnrm(lw) - pnorm_std(lw, 1L, 1L)),
              d_ub_ub = f_ub ? 0 : ub * d_ub,
              d_lb_lb = f_lb ? 0 : lw * d_lb,
               sd_inv = std::sqrt(*sig_inv);

    out[0L] = p_ub - p_lb;
    double const d_mu = -(d_ub - d_lb) * sd_inv;
    for(arma::uword j = 0; j < n_fix; ++j)
      out[j + 1] = X.at(0, j) * d_mu;

    double const d_sig = -(d_ub_ub - d_lb_lb) / 2 * sd_inv * sd_inv;
    for(unsigned s = 0; s  < scale_mats.size(); ++s)
      out[s + n_fix + 1] = d_sig * scale_mats.at(s).at(0, 0);
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
    size_t minvls;
    int inform;
    /// maximum estimated absolute accuracy of finest
    double abserr;
    /// likelihood approximation
    double likelihood;
    /// the derivative approximation
    arma::vec derivs;
    /// the approximate standard errors
    arma::vec sd_errs;
  };

  out_type get_output(double * res,  double const * sdest, size_t const minvls,
                      int const inform, double const abserr,
                      int const *indices){
    out_type out;
    out.minvls = minvls;
    out.inform = inform;
    out.abserr = abserr;

    double const likelihood = *res;
    out.likelihood = likelihood;
    out.sd_errs = arma::vec(sdest, get_n_integrands());

    // add terms to the derivative w.r.t. the scale parameters
    if(n_mem > 1){
      out.likelihood *= norm_const;
      out.sd_errs[0] *= norm_const;

      {
        // correct for the normalization constant
        double const rel_likelihood = norm_const / out.likelihood;
        for(unsigned i = 1; i <= n_fix + n_scales; ++i){
          res        [i] *= rel_likelihood;
          // we ignore the uncertainty from the likelihood approximation.
          // This is typically quite small compared to that of the derivative
          // of the likelihood
          out.sd_errs[i] *= rel_likelihood;
        }
      }

      double * PEDMOD_RESTRICT const d_sc = res + n_fix + 1L;
      double * sig_inv_ele = sig_inv;
      for(size_t s = 0; s < n_scales; ++s)
        scale_mats_ptr[s] = scale_mats.at(s).begin();

      for(arma::uword c = 0; c < n_mem; ++c){
        for(arma::uword r = 0; r < c; ++r, ++sig_inv_ele)
          for(size_t s = 0; s < n_scales; ++s)
            d_sc[s] -= *sig_inv_ele * *scale_mats_ptr[s]++;

        for(size_t s = 0; s < n_scales; ++s){
          d_sc[s] -= .5 * *sig_inv_ele * *scale_mats_ptr[s];
          scale_mats_ptr[s] += n_mem - c;
        }
        ++sig_inv_ele;
      }
    }

    // set deriv elements
    arma::vec &derivs = out.derivs;
    arma::uword const n_derivs = n_fix + n_scales;
    derivs.set_size(n_derivs);
    std::copy(res + 1, res + 1 + n_derivs, derivs.begin());
    return out;
  }
};

/**
 * Computes the derivatives for the integral
 *
 *   int_B phi(u; mu, Sigma) du
 *
 * The first element of the returned object is the likelihood factor. The
 * additional elements are the derivatives w.r.t. mu and Sigma with the latter
 * stored as the full k x k matrix ignoring the symmetry.
 */
class generic_l_factor {
  unsigned const n_vars,
           n_integrands = get_n_integrands(n_vars);

  /// working memory
  static cache_mem<double> dmem;

  double * cdf_mem() const { return dmem.get_mem(); }

  double * Sig_chol_tri() { return cdf_mem() + 2 * n_integrands; }

  double * internal_mem() {
    return Sig_chol_tri() + (n_vars * (n_vars + 1)) / 2;
  }

  /// the normalization constant
  double const norm_const;

  static unsigned get_n_integrands(unsigned const max_dim){
    return 1 + (max_dim * (max_dim + 3)) / 2;
  }

public:
  /// must be called prior to calling the member functions in the class
  static void alloc_mem
    (unsigned const max_dim, unsigned const max_threads,
     unsigned const max_n_sequences){
    rand_Korobov<cdf<generic_l_factor> >::alloc_mem(
        max_dim, get_n_integrands(max_dim), max_threads);
    sobol_wrapper<cdf<generic_l_factor> >::alloc_mem(
        max_dim, get_n_integrands(max_dim), max_threads, max_n_sequences);

    size_t n_mem = 2 * get_n_integrands(max_dim);
    n_mem += (max_dim * (max_dim + 1)) / 2;
    n_mem += max_dim * n_qmc_seqs();
    dmem.set_n_mem(n_mem, max_threads);
  }

  generic_l_factor(arma::vec const &mu, arma::mat const &Sig,
                   double const norm_const):
    n_vars{mu.n_elem}, norm_const{norm_const} {
      if(mu.n_elem != Sig.n_rows)
        throw std::invalid_argument("mu.n_elem != Sig.n_rows");
      else if(Sig.n_cols != Sig.n_rows)
        throw std::invalid_argument("Sig.n_cols != Sig.n_rows");
    }

  void prep_permutated(arma::mat const &Sig, int const *indices) {
    arma::mat Sig_chol = arma::chol(Sig);
    copy_upper_tri(Sig_chol, Sig_chol_tri());
  }

  unsigned get_n_integrands() PEDMOD_NOEXCEPT {
    return n_integrands;
  }

  double * get_wk_mem() PEDMOD_NOEXCEPT {
    return cdf_mem();
  }

  constexpr static bool needs_last_unif() PEDMOD_NOEXCEPT {
    return true;
  }

  double get_norm_constant() PEDMOD_NOEXCEPT {
    return norm_const;
  }

  void operator()
  (double const * PEDMOD_RESTRICT draw, double * PEDMOD_RESTRICT out,
   int const *, bool const, unsigned const n_draws){
    for(unsigned k = 0; k < n_draws; ++k)
      out[k * n_integrands] = 1;

    double * PEDMOD_RESTRICT const draw_scaled = internal_mem();
    for(unsigned v = 0; v < n_vars; ++v)
      for(unsigned k = 0; k < n_draws; ++k)
        draw_scaled[v + k * n_vars] = draw[k + v * n_draws];

    constexpr char uplo{'U'}, trans{'N'}, diag{'N'};
    constexpr int incx{1};
    int const int_n_vars = n_vars;
    for(unsigned k = 0; k < n_draws; ++k)
      F77_CALL(dtpsv)
        (&uplo, &trans, &diag, &int_n_vars, Sig_chol_tri(),
         draw_scaled + k * n_vars, &incx, 1, 1, 1);

    {
      double * PEDMOD_RESTRICT const d_mean{out + 1};
      for(unsigned k = 0; k < n_draws; ++k)
        std::copy(draw_scaled + k * n_vars, draw_scaled + (k + 1) * n_vars,
                  d_mean + k * n_integrands);
    }

    double * PEDMOD_RESTRICT const d_vcov{out + 1 + n_vars};
    for(unsigned k = 0; k < n_draws; ++k){
      size_t upper_idx{};
      for(unsigned v1 = 0; v1 < n_vars; ++v1)
        for(unsigned v2 = 0; v2 <= v1; ++v2, ++upper_idx)
          d_vcov[upper_idx + k * n_integrands] =
            draw_scaled[v1 + k * n_vars] * draw_scaled[v2 + k * n_vars];
    }
  }

  void univariate(double * out, double const lw, double const ub){
    constexpr double log_sqrt_2_pi_inv{0.918938533204673};
    auto log_dnrm = [&](double const x){
      return -x * x / 2. - log_sqrt_2_pi_inv;
    };

    bool const f_ub = std::isinf(ub),
               f_lb = std::isinf(lw);

    double const p_ub = f_ub ? 1 : pnorm_std(ub, 1L, 0L),
                 p_lb = f_lb ? 0 : pnorm_std(lw, 1L, 0L),
                 d_ub = f_ub ? 0 : std::exp(log_dnrm(ub) - pnorm_std(ub, 1L, 1L)),
                 d_lb = f_lb ? 0 : std::exp(log_dnrm(lw) - pnorm_std(lw, 1L, 1L)),
                 d_ub_ub = f_ub ? 0 : ub * d_ub,
                 d_lb_lb = f_lb ? 0 : lw * d_lb,
                 sd_inv = 1 / *Sig_chol_tri();

    out[0] = p_ub - p_lb;
    out[1] = -(d_ub - d_lb) * sd_inv;
    out[2] = -(d_ub_ub - d_lb_lb) / 2 * sd_inv * sd_inv;
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
    size_t minvls;
    int inform;
    /// maximum estimated absolute accuracy of finest
    double abserr;
    /// likelihood approximation
    double likelihood;
    /// the derivative approximation
    arma::vec derivs;
    /// the approximate standard errors
    arma::vec sd_errs;
  };

  out_type get_output(double * res,  double const * sdest, size_t const minvls,
                      int const inform, double const abserr,
                      int const *indices){
    out_type out;
    out.minvls = minvls;
    out.inform = inform;
    out.abserr = abserr;

    double const likelihood = *res;
    out.likelihood = likelihood;
    out.sd_errs = arma::vec(sdest, get_n_integrands());

    if(n_vars > 1){
      // setup the derivative w.r.t. the covariance matrix and account for the
      // permutation

      out.likelihood *= norm_const;
      out.sd_errs[0] *= norm_const;

      {
        // correct for the normalization constant and get the derivative w.r.t.
        // the log likelihood
        double const rel_likelihood = norm_const / out.likelihood;
        std::for_each(res + 1, res + n_integrands,
                      [&](double &x){ x *= rel_likelihood; });
        std::for_each(out.sd_errs.begin() + 1, out.sd_errs.end(),
                      [&](double &x){ x *= rel_likelihood; });
      }

      out.derivs.resize(n_vars * (n_vars + 1));

      arma::vec d_mean(out.derivs.begin(), n_vars, false);

      {
        double const * const d_mean_permu{res + 1};
        for(unsigned v = 0; v < n_vars; ++v)
          d_mean[indices[v]] = d_mean_permu[v];
      }

      arma::mat d_Sig(d_mean.end(), n_vars, n_vars, false);

      std::unique_ptr<double[]> Sig_permu_inv
        (new double[(n_vars * (n_vars + 1)) / 2]);

      {
        std::copy
          (Sig_chol_tri(), Sig_chol_tri() + (n_vars * (n_vars + 1)) / 2,
           Sig_permu_inv.get());
        constexpr char uplo{'U'};
        int const int_n_vars = n_vars;
        int info{};
        F77_CALL(dpptri)(&uplo, &int_n_vars, Sig_permu_inv.get(), &info, 1);

        if(info != 0)
          throw std::runtime_error("dpptri failed");
      }

      double const * const d_Sig_permu_outer{res + 1 + n_vars};
      size_t idx_upper{};
      for(unsigned v1 = 0; v1 < n_vars; ++v1)
        for(unsigned v2 = 0; v2 <= v1; ++v2, ++idx_upper){
          double const deriv_val
           {(d_Sig_permu_outer[idx_upper] - Sig_permu_inv[idx_upper]) / 2};

          d_Sig(indices[v2], indices[v1]) = deriv_val;
          d_Sig(indices[v1], indices[v2]) = deriv_val;
        }

      return out;
    }

    // set deriv elements
    arma::vec &derivs = out.derivs;
    arma::uword const n_derivs = n_integrands - 1;
    derivs.set_size(n_derivs);
    std::copy(res + 1, res + 1 + n_derivs, derivs.begin());
    return out;
  }
};

} // namespace pedmod

#endif
