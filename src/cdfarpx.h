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
 * @param X Matrix to copy.
 * @param x Pointer to copy to.
 */
inline void copy_upper_tri
  (arma::mat const &X, double * __restrict__ x) noexcept {
  int const p = X.n_cols;
  for(int c = 0; c < p; c++)
    for(int r = 0; r <= c; r++, x++)
      *x = X.at(r, c);
}

/**
 * copies the lower triangular matrix.
 *
 * @param X Matrix to copy.
 * @param x Pointer to copy to.
 */
inline void copy_lower_tri
  (arma::mat const &X, double * __restrict__ x) noexcept {
  int const p = X.n_cols;
  for(int c = 0; c < p; c++)
    for(int r = c; r < p; r++, x++)
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
  const int ndim,
            n_integrands;
  const bool use_aprx;
  bool is_permutated = false;

  static constexpr bool const
    needs_last_unif = T_Functor::needs_last_unif();

  // cached memory to use
  static cache_mem<int   > imem;
  static cache_mem<double> dmem;

  arma::ivec infin;
  arma::ivec indices;

  double * __restrict__ const lower      = dmem.get_mem(),
         * __restrict__ const upper      = lower + ndim,
         * __restrict__ const sigma_chol = upper + ndim,
         * __restrict__ const draw       =
         sigma_chol + (ndim * (ndim + 1L)) / 2L,
         * __restrict__ const dtmp_mem   = draw + ndim;

  // memory that can be used
  int * const itmp_mem = indices.end();

public:
  /**
   * must be called prior to calling the constructor or any member
   * functions.
   */
  static void alloc_mem(int const max_ndim, int const max_threads) {
    int const n_up_tri = (max_ndim * (max_ndim + 1)) / 2;

    imem.set_n_mem(3 * max_ndim                                 ,
                   max_threads);
    dmem.set_n_mem(7 * max_ndim + n_up_tri + max_ndim * max_ndim,
                   max_threads);
  }

  cdf(T_Functor &functor, arma::vec const &lower_in,
      arma::vec const &upper_in, arma::vec const &mu_in,
      arma::mat const &sigma_in, bool const do_reorder,
      bool const use_aprx):
    functor(functor),
    ndim(mu_in.n_elem),
    n_integrands(functor.get_n_integrands()),
    use_aprx(use_aprx),
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
    auto get_dmem = [&](int const n_ele) -> double * {
      double * out = cur_dtmp_mem;
      cur_dtmp_mem += n_ele;
      return out;
    };

    double * sds = get_dmem(ndim);
    for(int i = 0; i < ndim; ++i){
      sds[i] = std::sqrt(sigma_in.at(i, i));
      lower[i] = (lower_in[i] - mu_in[i]) / sds[i];
      upper[i] = (upper_in[i] - mu_in[i]) / sds[i];
    }

    is_permutated = false;
    {
      int *idx = indices.begin();
      for(int i = 0; i < ndim; ++i, ++idx)
        *idx = i;
    }

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
             nddim = ndim;
      std::fill(delta, delta + ndim, 0.);
      arma::ivec infi(itmp_mem, ndim, false);

      F77_CALL(mvsort)(
        &ndim, lower, upper, delta,
        correl.cor_vec.memptr(), infin.begin(), y, &pivot, &nddim, A, B,
        DL, sigma_chol, infi.memptr(), &F_inform, indices.begin(),
        &doscale);

      if(F_inform != 0)
        throw std::runtime_error("cdf::cdf: error in mvsort");

      for(int i = 0; i < ndim; ++i){
        if(indices[i] != i){
          is_permutated = true;
          break;
        }
      }

      if(is_permutated){
        for(int i = 0; i < ndim; ++i){
          lower[i] = A[i];
          upper[i] = B[i];
          infin[i] = infi[i];
        }

        arma::mat sigma_permu(get_dmem(ndim * ndim), ndim, ndim, false);
        for(int j = 0; j < ndim; ++j)
          for(int i = 0; i < ndim; ++i)
            sigma_permu.at(i, j) = sigma_in.at(indices[i], indices[j]);
        functor.prep_permutated(sigma_permu);

        return;

      } else
        for(int i = 0; i < ndim; ++i){
          lower[i] = A[i];
          upper[i] = B[i];
        }

    } else if(!do_reorder and ndim > 1L) {
      arma::mat tmp(get_dmem(ndim * ndim), ndim, ndim, false);
      tmp = sigma_in;
      for(int i = 0; i < ndim; ++i)
        for(int j = 0; j <ndim; ++j)
          tmp.at(i, j) /= sds[i] * sds[j];

      if(!arma::chol(tmp, tmp)) // TODO: memory allocation
        std::fill(sigma_chol, sigma_chol + (ndim * (ndim + 1L)) / 2L,
                  std::numeric_limits<double>::infinity());
      else
        copy_upper_tri(tmp, sigma_chol);

      if(ndim > 1L){
        /* rescale such that Choleksy decomposition has ones in the diagonal */
        double * sc = sigma_chol;
        for(int i = 0; i < ndim; ++i){
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
           * const __restrict__ dr  = draw;

    double w(1.);
    double const * __restrict__ sc   = sigma_chol,
                 * __restrict__ lw   = lower,
                 * __restrict__ up   = upper,
                 * __restrict__ unif = unifs;
    int const *infin_j = infin.begin();
    /* loop over variables and transform them to truncated normal
     * variables */
    for(int j = 0; j < ndim; ++j, ++sc, ++lw, ++up, ++infin_j, ++unif){
      double su(0.);
      for(int i = 0; i < j; ++i, sc++)
        su += *sc * dr[i];

      auto pnorm_use = [&](double const x) -> double {
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
        double const l_diff = lim_u - lim_l;
        w *= l_diff;

        if(needs_last_unif or j + 1 < ndim){
          double const quant_val = lim_l + *unif * l_diff;
          dr[j] =
            use_aprx ?
            safe_qnorm_aprx(quant_val) :
            safe_qnorm_w   (quant_val);
        }

      } else {
        w = 0;
        std::fill(dr + j, dr + ndim, 0.);
        break;

      }
    }

    /* evaluate the integrand and weight the result. */
    functor(dr, out, indices.begin(), is_permutated);

    w /= functor.get_norm_constant();
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
   cdf_methods const method, int const minvls, unsigned const n_sequences){
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

  constexpr static int get_n_integrands() {
    return 1L;
  }

  inline void operator()
    (double const *, double * out, int const *, bool const)
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
    int minvls, inform;
    /// maximum norm of estimated absolute accuracy of finest
    double abserr;
    /// likelihood approximation
    double likelihood;
  };

  out_type get_output(double const * res, double const * sdest,
                      int const minvls, int const inform, double const abserr,
                      int const *){
    return out_type { minvls, inform, abserr, *res };
  }

  inline void prep_permutated(arma::mat const&) { }
};

/**
 * functor classes used as template argument for cdf used to approximate the
 * derivatives of the likelihood factors for each family.
 *
 * The returned approximations is a) the likelihood factor, b) the
 * derivative of log likelihood w.r.t. the mean vector, and w.r.t. each of the
 * scale parameters.
 */
class pedigree_l_factor {
public:
  /// the scale matrices for the different effects
  std::vector<arma::mat> const scale_mats;
  /// the number of members in this family
  int const n_mem = scale_mats[0].n_rows;
  /// design matrix for the fixed effects
  arma::mat const X;
  /// the number of fixed effects
  int const n_fix = X.n_cols,
     n_integrands = 1 + n_fix + scale_mats.size(),
         n_scales = scale_mats.size();
  /// should sparse matrices be used
  bool const use_sparse;
  /// the scale matrices as sparse matrices
  std::vector<arma::sp_mat> scale_mats_sparse =
    ([&]() -> std::vector<arma::sp_mat> {
      std::vector<arma::sp_mat> out;
      if(use_sparse){
        out.reserve(scale_mats.size());
        for(arma::mat const &m : scale_mats){
          arma::mat tmp = arma::trimatu(m);
          tmp.diag() *= .5;
          out.emplace_back(std::move(tmp));
        }
      }

      return out;
    })();

private:
  /// working memory
  static cache_mem<double> dmem;

  /**
   * points to the upper triangular part of the inverse of the Cholesky
   * decomposition.
   */
  double * sigma_chol_inv;

  /// points to the upper triangular part of the inverse.
  double * sig_inv;

  /// working memory to be used by cdf
  double * cdf_mem;

  /// array of pointer to the scale matrices' element which we will need.
  std::unique_ptr<double const *[]> scale_mats_ptr =
    std::unique_ptr<double const *[]>(new double const *[n_scales]);

  /// the normalization constant
  double norm_const = std::numeric_limits<double>::quiet_NaN();

public:
  /// sets the scale matrices. There are no checks on the validity
  pedigree_l_factor(std::vector<arma::mat> const &scale_mats,
                    unsigned const max_threads, arma::mat const &X_in,
                    bool const use_sparse, unsigned const max_n_sequences):
  scale_mats(scale_mats), X(X_in.t()), use_sparse(use_sparse) {
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
    dmem.set_n_mem(
      2 * n_mem * n_mem + n_mem * (n_mem + 1) + 2 * get_n_integrands() +
        2 * n_mem,
      max_threads);
  }

  inline int get_n_integrands() PEDMOD_NOEXCEPT {
    return n_integrands;
  }

  inline double * get_wk_mem() PEDMOD_NOEXCEPT {
    return cdf_mem;
  }

  constexpr static bool needs_last_unif() PEDMOD_NOEXCEPT {
    return true;
  }

  inline double get_norm_constant() PEDMOD_NOEXCEPT {
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

    // create the objects we need
    arma::mat t1(dmem.get_mem(), n_mem, n_mem, false),
              t2(t1.end()      , n_mem, n_mem, false);
    if(!arma::chol(t1, sig))
      throw std::runtime_error("pedigree_ll_factor::setup: chol failed");
    if(!arma::inv(t2, t1))
      throw std::runtime_error("pedigree_ll_factor::setup: inv failed");
    sigma_chol_inv = t2.end();
    copy_upper_tri(t2, sigma_chol_inv);

    if(!arma::inv_sympd(t1, sig))
      throw std::runtime_error("pedigree_ll_factor::setup: inv_sympd failed");
    sig_inv = sigma_chol_inv + (n_mem * (n_mem + 1)) / 2;
    copy_upper_tri(t1, sig_inv);

    cdf_mem = sig_inv + (n_mem * (n_mem + 1)) / 2;
  }

  void prep_permutated(arma::mat const &sigma_permu) {
    // need to re-compute the inverse of the Cholesky decomposition
    arma::mat t1(dmem.get_mem(), n_mem, n_mem, false),
              t2(t1.end()      , n_mem, n_mem, false);
    if(!arma::chol(t1, sigma_permu))
      throw std::runtime_error("pedigree_ll_factor::setup: chol failed");
    if(!arma::inv(t2, t1))
      throw std::runtime_error("pedigree_ll_factor::setup: inv failed");
    copy_upper_tri(t2, sigma_chol_inv);
  }

  inline void operator()
    (double const * __restrict__ draw, double * __restrict__ out,
     int const *indices, bool const is_permutated) {
      *out = 1;
      double * __restrict__ const d_mu      = get_wk_mem() + 2 * get_n_integrands(),
             * __restrict__ const d_mu_perm = d_mu + n_mem,
             * __restrict__ const d_fix     = out + 1,
             * __restrict__ const d_sc      = d_fix + n_fix;

      std::fill(d_mu   , d_mu + n_mem            , 0.);
      std::fill(out + 1, out + get_n_integrands(), 0.);

      // handle the derivatives w.r.t. the mean
      {
        double *sig_chov_inv_ele = sigma_chol_inv;
        double const * d = draw;
        for(int c = 0; c < n_mem; ++c, ++d){
          double *rhs = d_mu;
          for(int r = 0; r <= c; ++r, ++rhs, ++sig_chov_inv_ele)
            *rhs += *d * *sig_chov_inv_ele;
        }
      }

      // permute back
      for(int i = 0; i < n_mem; ++i)
        d_mu_perm[indices[i]] = d_mu[i];

      // derivatives w.r.t. the fixed effects
      {
        double const *xij = X.begin();
        for(int j = 0; j < n_fix; ++j, xij += n_mem){
          for(int i = 0; i < n_mem; ++i)
            d_fix[j] += xij[i] * d_mu_perm[i];
        }
      }

      // handle the derivatives w.r.t. the scale parameters
      if(use_sparse){
        for(int s = 0; s < n_scales; ++s){
          arma::sp_mat::const_iterator it = scale_mats_sparse[s].begin();
          arma::sp_mat::const_iterator const end = scale_mats_sparse[s].end();
          for(; it != end; ++it)
            d_sc[s] += d_mu_perm[it.col()] * d_mu_perm[it.row()] * *it;
        }

        return;
      }

      for(int s = 0; s < n_scales; ++s)
        scale_mats_ptr[s] = scale_mats.at(s).begin();

      for(int c = 0; c < n_mem; ++c){
        double * __restrict__ d_mu_prod = d_mu;
        for(int r = 0; r < c; ++r)
          d_mu_prod[r] = d_mu_perm[r] * d_mu_perm[c];
        d_mu_prod[c] = .5 * d_mu_perm[c] * d_mu_perm[c];

        for(int s = 0; s < n_scales; ++s)
          for(int r = 0; r <= c; ++r)
            d_sc[s] += d_mu_prod[r] * scale_mats_ptr[s][r];

        for(int s = 0; s < n_scales; ++s)
          scale_mats_ptr[s] += n_mem;
      }
    }

  inline void univariate(double * out, double const lw, double const ub) {
    constexpr double const log_sqrt_2_pi_inv = 0.918938533204673;
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
               sd_inv = *sigma_chol_inv;

    out[0L] = p_ub - p_lb;
    double const d_mu = -(d_ub - d_lb) * sd_inv;
    for(int j = 0; j < n_fix; ++j)
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
    int minvls, inform;
    /// maximum estimated absolute accuracy of finest
    double abserr;
    /// likelihood approximation
    double likelihood;
    /// the derivative approximation
    arma::vec derivs;
    /// the approximate standard errors
    arma::vec sd_errs;
  };

  out_type get_output(double * res,  double const * sdest, int const minvls,
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
        for(int i = 1; i <= n_fix + n_scales; ++i){
          res        [i] *= rel_likelihood;
          // we ignore the uncertainty from the likelihood approximation.
          // This is typically quite small compared to that of the derivative
          // of the likelihood
          out.sd_errs[i] *= rel_likelihood;
        }
      }

      double * __restrict__ const d_sc = res + n_fix + 1L;
      double * sig_inv_ele = sig_inv;
      for(int s = 0; s < n_scales; ++s)
        scale_mats_ptr[s] = scale_mats.at(s).begin();

      for(int c = 0; c < n_mem; ++c){
        for(int r = 0; r < c; ++r, ++sig_inv_ele)
          for(int s = 0; s < n_scales; ++s)
            d_sc[s] -= *sig_inv_ele * *scale_mats_ptr[s]++;

        for(int s = 0; s < n_scales; ++s){
          d_sc[s] -= .5 * *sig_inv_ele * *scale_mats_ptr[s];
          scale_mats_ptr[s] += n_mem - c;
        }
        ++sig_inv_ele;
      }
    }

    // set deriv elements
    arma::vec &derivs = out.derivs;
    int const n_derivs = n_fix + n_scales;
    derivs.set_size(n_derivs);
    std::copy(res + 1, res + 1 + n_derivs, derivs.begin());
    return out;
  }
};
} // namespace pedmod

#endif
