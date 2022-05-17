#include "cdfaprx.h"
#include <R_ext/RS.h>

namespace pedmod {

extern "C" {
  /// either performs a forward or backward solve
  void F77_NAME(dtpsv)
    (const char * /* uplo */, const char * /* trans */, const char * /* diag */,
     const int * /* n */, const double * /* ap */, double * /* x */,
     const int* /* incx */, size_t, size_t, size_t);
}

arma::ivec get_infin(
    arma::ivec &out, arma::vec const &lower, arma::vec const &upper){
#ifdef DO_CHECKS
  arma::uword const n = lower.size();
  if(upper.size() != n)
    throw std::invalid_argument("get_infin: invalid 'upper'");
#endif
  double const *l = lower.begin(),
               *u = upper.begin();
  for(auto &o : out){
    bool const li = std::isinf(*l++),
               ui = std::isinf(*u++);
    if      ( li and  ui)
      o = -1L;
    else if ( li and !ui)
      o =  0L;
    else if (!li and  ui)
      o =  1L;
    else
      o =  2L;
  }

  return out;
}

cor_vec_res get_cor_vec(const arma::mat &cov){
  cor_vec_res out;
  arma::vec     &sds = out.sds,
            &cor_vec = out.cor_vec;

  arma::uword const n = cov.n_cols;
  sds = arma::sqrt(cov.diag());
  cor_vec.resize((n * (n - 1L)) / 2L);

#ifdef DO_CHECKS
  if(n != cov.n_rows)
    throw std::invalid_argument("get_cor_vec: invalid 'cov'");
  if(n <= 0L)
    throw std::invalid_argument("get_cor_vec: invalid 'n'");
#endif

  double *o = cor_vec.begin();
  for(unsigned c = 1L; c < n; ++c)
    for(unsigned r = 0; r < c; ++r)
      *o++ = cov(r, c) / sds[c] / sds[r];

  return out;
}

pedigree_l_factor::pedigree_l_factor
  (std::vector<arma::mat> const &scale_mats, unsigned const max_threads,
   arma::mat const &X_in, unsigned const max_n_sequences):
  scale_mats(scale_mats), X(X_in.t()) {
  // checks
  if(scale_mats.size() < 1)
    throw std::invalid_argument("pedigree_l_factor::pedigree_l_factor: no scale matrices are passed");
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

void pedigree_l_factor::setup
  (arma::mat &sig, double const *scales, double const norm_constant_arg,
   bool const only_cov){
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

void pedigree_l_factor::prep_permutated
  (arma::mat const &sig, int const *indices) {
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

void pedigree_l_factor::operator()
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

void pedigree_l_factor::univariate
  (double * out, double const lw, double const ub) {
  constexpr double log_sqrt_2_pi_inv{0.918938533204673};
  auto log_dnrm = [&](double const x){
    return -x * x / 2. - log_sqrt_2_pi_inv;
  };

  // TODO: the code does not seem correct if one of the limits are not Inf
  bool const f_ub = std::isinf(ub),
             f_lb = std::isinf(lw);

  double const p_ub = f_ub ? 1 : pnorm_std(ub, 1L, 0L),
               p_lb = f_lb ? 0 : pnorm_std(lw, 1L, 0L),
               d_ub = f_ub ? 0 : std::exp(log_dnrm(ub) - pnorm_std(ub, 1L, 1L)),
               d_lb = f_lb ? 0 : std::exp(log_dnrm(lw) - pnorm_std(-lw, 1L, 1L)),
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

pedigree_l_factor::out_type pedigree_l_factor::get_output
  (double * res,  double const * sdest, size_t const minvls,
   int const inform, double const abserr, int const *indices){
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

pedigree_l_factor_Hessian::pedigree_l_factor_Hessian
  (std::vector<arma::mat> const &scale_mats, unsigned const max_threads,
   arma::mat const &X_in, unsigned const max_n_sequences):
  scale_mats(scale_mats), X(X_in.t()) {
  // checks
  if(scale_mats.size() < 1)
    throw std::invalid_argument("pedigree_l_factor_Hessian::pedigree_l_factor_Hessian: no scale matrices are passed");
  arma::uword const u_mem = n_mem;
  for(auto &S : scale_mats)
    if(S.n_rows != u_mem or S.n_rows != u_mem)
      throw std::invalid_argument("pedigree_l_factor_Hessian::pedigree_l_factor_Hessian: not all scale matrices are square matrices or have the same dimensions");
  if(X.n_rows != u_mem)
    throw std::invalid_argument("pedigree_l_factor_Hessian::pedigree_l_factor_Hessian: invalid X");

  // setup working memory
  rand_Korobov<cdf<pedigree_l_factor_Hessian> >::alloc_mem(
      n_mem, get_n_integrands(), max_threads);
  sobol_wrapper<cdf<pedigree_l_factor_Hessian> >::alloc_mem(
      n_mem, get_n_integrands(), max_threads, max_n_sequences);

  size_t const working_memory{n_fix + n_scales};
  dmem.set_n_mem(
    (n_mem * (n_mem + 1)) / 2 + // vcov_chol
      n_mem * n_mem + // vcov_inv
      n_mem * n_fix + // X_permu
      n_mem * n_mem * n_scales + // scale_mats_permu
      2 * get_n_integrands() +
      working_memory,
      max_threads);
}

void pedigree_l_factor_Hessian::setup
  (arma::mat &sig, double const *scales, double const norm_constant_arg){
  sig.zeros(n_mem, n_mem);
  sig.diag() += 1;
  for(unsigned i = 0; i < scale_mats.size(); ++i)
    sig += scales[i] * scale_mats[i];
  norm_const = norm_constant_arg;

  vcov_chol = dmem.get_mem();
  vcov_inv = vcov_chol + (n_mem * (n_mem + 1)) / 2;
  X_permu = vcov_inv + n_mem * n_mem;

  double * next_mem = X_permu + n_mem * n_fix;
  for(arma::uword i = 0; i < n_scales; ++i, next_mem += n_mem * n_mem)
    scale_mats_permu[i] = next_mem;

  cdf_mem = next_mem;
  interal_mem = cdf_mem + 2 * get_n_integrands();
}

void pedigree_l_factor_Hessian::prep_permutated
  (arma::mat const &sig, int const *indices) {
  arma::mat const sig_chol = arma::chol(sig);

  {
    double * vcov_chol_id{vcov_chol};
    for(arma::uword id = 0; id < n_mem; ++id, vcov_chol_id += id)
      std::copy(sig_chol.colptr(id), sig_chol.colptr(id) + id + 1, vcov_chol_id);
  }

  {
    arma::mat const sig_inv = arma::inv_sympd(sig);
    std::copy(sig_inv.begin(), sig_inv.end(), vcov_inv);
  }

  for(arma::uword fix = 0; fix < n_fix; ++fix)
    for(arma::uword id = 0; id < n_mem; ++id)
      X_permu[id + fix * n_mem] = X(indices[id], fix);

  for(size_t scale = 0; scale < n_scales; ++scale)
    for(arma::uword id1 = 0; id1 < n_mem; ++id1)
      for(arma::uword id2 = 0; id2 < n_mem; ++id2)
        *(scale_mats_permu[scale] + id2 + id1 * n_mem) =
          scale_mats[scale](indices[id2], indices[id1]);
}

void pedigree_l_factor_Hessian::operator()
  (double const * PEDMOD_RESTRICT draw, double * PEDMOD_RESTRICT out,
   int const *, bool const, unsigned const n_draws) {
  for(unsigned idx_draw = 0; idx_draw < n_draws; ++idx_draw)
    out[idx_draw * n_integrands] = 1;

  size_t const shift_scaled_res{1},
               shift_outer_res{shift_scaled_res + n_mem},
               shift_vec_outer_res{shift_outer_res + n_mem * n_mem};

  /*
   * let
   *
   *   Psi be the covariance matrix
   *   S^TS = Psi be the Cholesky decomposition
   *   U be the draws
   *
   * then we first compute S^(-1)u and (S^(-1)uu^TS^(-1) - Psi^(-1)) / 2
   */

  double * const outer_vec{interal_mem};
  for(unsigned idx_draw = 0; idx_draw < n_draws; ++idx_draw){
    // fill in S^(-1)u and (S^(-1)uu^TS^(-1) - Psi^(-1)) / 2
    double * const res_i{out + idx_draw * get_n_integrands()};
    double * const draw_scaled{res_i + shift_scaled_res};
    {
      for(unsigned id = 0; id < n_mem; ++id)
        draw_scaled[id] = draw[id * n_draws + idx_draw];
      constexpr char uplo{'U'},
                    trans{'N'},
                     diag{'N'};
      int const n_mem_i = n_mem;
      constexpr int incx{1};
      F77_CALL(dtpsv)
        (&uplo, &trans, &diag, &n_mem_i, vcov_chol, draw_scaled,
         &incx, 1, 1, 1);
    }

    double * const outer_res{res_i + shift_outer_res};
    for(unsigned id1 = 0; id1 < n_mem; ++id1)
      for(unsigned id2 = 0; id2 < n_mem; ++id2)
        outer_res[id2 + id1 * n_mem] =
          (draw_scaled[id1] * draw_scaled[id2]
             - vcov_inv[id2 + id1 * n_mem]) / 2;

    /// compute the required outer product for the hessian
    for(unsigned fix = 0; fix < n_fix; ++fix)
      outer_vec[fix] = std::inner_product
        (X_permu + fix * n_mem, X_permu + (fix + 1) * n_mem, draw_scaled, 0.);

    for(size_t scale = 0; scale < n_scales; ++scale)
      outer_vec[scale + n_fix] = std::inner_product
        (outer_res, outer_res + n_mem * n_mem, scale_mats_permu[scale], 0.);

    size_t const vec_outer_res_dim{n_fix + n_scales};
    double * const vec_outer_res{res_i + shift_vec_outer_res};
    for(size_t param1 = 0; param1 < vec_outer_res_dim; ++param1)
      for(size_t param2 = 0; param2 < vec_outer_res_dim; ++param2)
        vec_outer_res[param2 + param1 * vec_outer_res_dim] =
          outer_vec[param1] * outer_vec[param2];
  }
}

void pedigree_l_factor_Hessian::univariate
  (double * out, double const lw, double const ub) {
  constexpr double log_sqrt_2_pi_inv{0.918938533204673};
  auto log_dnrm = [&](double const x){
    return -x * x / 2. - log_sqrt_2_pi_inv;
  };

  // TODO: the code does not seem correct if one of the limits are not Inf
  bool const f_ub = std::isinf(ub),
             f_lb = std::isinf(lw);

  double const p_ub = f_ub ? 1 : pnorm_std(ub, 1L, 0L),
               p_lb = f_lb ? 0 : pnorm_std(lw, 1L, 0L),
               d_ub = f_ub ? 0 : std::exp(log_dnrm(ub) - pnorm_std(ub, 1L, 1L)),
               d_lb = f_lb ? 0 : std::exp(log_dnrm(lw) - pnorm_std(-lw, 1L, 1L)),
            d_ub_ub = f_ub ? 0 : ub * d_ub,
            d_lb_lb = f_lb ? 0 : lw * d_lb,
             sd_inv = std::sqrt(*vcov_inv);

  out[0L] = p_ub - p_lb;

  // the gradient
  double const d_mu = -(d_ub - d_lb) * sd_inv;
  for(arma::uword j = 0; j < n_fix; ++j)
    out[j + 1] = X.at(0, j) * d_mu;

  double const d_sig = -(d_ub_ub - d_lb_lb) / 2 * sd_inv * sd_inv;
  for(size_t s = 0; s  < n_scales; ++s)
    out[s + n_fix + 1] = d_sig * scale_mats.at(s).at(0, 0);

  arma::uword const hess_dim{n_fix + n_scales};
  double * const hess_begin{out + 1 + hess_dim};
  std::fill(hess_begin, hess_begin + hess_dim * hess_dim, 0);


  auto add_to_hess = [&](double const limit, double const ratio,
                         bool const is_upper){
    double const
        d_num
        {is_upper
          ? -(limit * ratio + ratio * ratio) * *vcov_inv
          : (limit * ratio - ratio * ratio) * *vcov_inv},
        d_cross
        {is_upper
          ? (-limit * ratio * ratio + ratio - limit * limit * ratio) *
              *vcov_inv * sd_inv / 2
          : (-limit * ratio * ratio - ratio + limit * limit * ratio) *
              *vcov_inv * sd_inv / 2},
        d_denom
        {is_upper
          ? (-limit * limit * limit * ratio - limit * limit * ratio * ratio +
              3 * limit * ratio) * *vcov_inv * *vcov_inv / 4
          : (limit * limit * limit * ratio - limit * limit * ratio * ratio -
            3 * limit * ratio) * *vcov_inv * *vcov_inv / 4};

    for(arma::uword j = 0; j < n_fix; ++j)
      for(arma::uword i = 0; i < n_fix; ++i)
        hess_begin[i + j * hess_dim] += d_num * X.at(0, j) * X.at(0, i);

    for(arma::uword s = 0; s < n_scales; ++s)
      for(arma::uword i = 0; i < n_fix; ++i)
        hess_begin[i + (s + n_fix) * hess_dim] +=
          d_cross * X.at(0, i) * scale_mats.at(s).at(0, 0);

    for(arma::uword s = 0; s < n_scales; ++s)
      for(arma::uword k = 0; k < n_scales; ++k)
        hess_begin[k + n_fix + (s + n_fix) * hess_dim] +=
          d_denom * scale_mats.at(k).at(0, 0) * scale_mats.at(s).at(0, 0);
   };

  if(!f_ub)
    add_to_hess(ub, d_ub, true);
  if(!f_lb)
    add_to_hess(lw, d_lb, false);

   arma::mat hess(hess_begin, hess_dim, hess_dim, false);
   hess = arma::symmatu(hess);
}

pedigree_l_factor_Hessian::out_type pedigree_l_factor_Hessian::get_output
  (double * res,  double const * sdest, size_t const minvls,
   int const inform, double const abserr, int const *indices){
  out_type out;
  out.minvls = minvls;
  out.inform = inform;
  out.abserr = abserr;

  size_t const hess_dim{n_fix + n_scales};
  arma::vec &gr = out.gradient;
  arma::mat &hess = out.hessian;
  out.sd_errs.resize(1 + hess_dim * (hess_dim + 1));

  if(n_mem > 1){
    size_t const shift_scaled_res{1},
                 shift_outer_res{shift_scaled_res + n_mem},
                 shift_vec_outer_res{shift_outer_res + n_mem * n_mem};

    out.sd_errs[0] = sdest[0] * norm_const;
    std::fill
      (out.sd_errs.begin() + 1, out.sd_errs.end(),
       std::numeric_limits<double>::quiet_NaN());
    out.likelihood = *res * norm_const;

    double const rel_likelihood = norm_const / out.likelihood;

    // compute the gradient
    gr.zeros(hess_dim);

    for(arma::uword fix = 0; fix < n_fix; ++fix)
      gr[fix] += std::inner_product
        (X_permu + fix * n_mem, X_permu + (fix + 1) * n_mem,
         res + shift_scaled_res, 0.);

    for(size_t scale = 0; scale < n_scales; ++scale)
      gr[scale + n_fix] += std::inner_product
        (res + shift_outer_res, res + shift_outer_res + n_mem * n_mem,
         scale_mats_permu[scale], 0.);

    // compute the Hessian
    arma::mat X_permu_mat(X_permu, n_mem, n_fix, false),
             vcov_inv_mat(vcov_inv, n_mem, n_mem, false); // TODO: solve instead

    hess.resize(hess_dim, hess_dim);
    std::copy
      (res + shift_vec_outer_res,
       res + shift_vec_outer_res + hess_dim * hess_dim, hess.begin());

    std::vector<arma::mat> scale_mats_arma;
    scale_mats_arma.reserve(n_scales);
    for(auto permu_mat : scale_mats_permu)
      scale_mats_arma.emplace_back(permu_mat, n_mem, n_mem, false);

    hess.submat(0, 0, n_fix - 1, n_fix - 1) -=
      (X_permu_mat.t() * vcov_inv_mat * X_permu_mat) / rel_likelihood;

    {
      arma::vec cross_vec(res + shift_scaled_res, n_mem, false);
      for(size_t scale = 0; scale < n_scales; ++scale)
        hess.submat(0, n_fix + scale, n_fix - 1, n_fix + scale) -=
          X_permu_mat.t() * vcov_inv_mat * scale_mats_arma[scale] * cross_vec;
    }
    {
      arma::mat outer_vec(res + shift_outer_res, n_mem, n_mem, false);
      arma::mat prod1;

      for(size_t scale1 = 0; scale1 < n_scales; ++scale1){
        prod1.zeros(n_mem, n_mem);
        prod1 += vcov_inv_mat * scale_mats_arma[scale1] * outer_vec;
        prod1 += outer_vec * scale_mats_arma[scale1] * vcov_inv_mat;
        prod1 +=
          vcov_inv_mat * scale_mats_arma[scale1] * vcov_inv_mat /
          (2 * rel_likelihood);

        for(size_t scale2 = 0; scale2 <= scale1; ++scale2)
          hess(n_fix + scale2, n_fix + scale1) -=
            std::inner_product
              (prod1.begin(), prod1.end(), scale_mats_permu[scale2], 0.);
      }
    }

    gr *= rel_likelihood;
    hess *= rel_likelihood;
    hess -= gr * gr.t();
    hess = arma::symmatu(hess);

    return out;
  }

  out.likelihood = *res;
  gr.resize(hess_dim);
  hess.resize(hess_dim, hess_dim);
  std::copy(res + 1, res + 1 + hess_dim, gr.begin());
  std::copy
    (res + 1 + hess_dim, res + 1 + hess_dim * (1 + hess_dim), hess.begin());
  out.sd_errs.zeros();

  return out;
}

void generic_l_factor::alloc_mem
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

void generic_l_factor::operator()
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

void generic_l_factor::univariate(double * out, double const lw, double const ub){
  constexpr double log_sqrt_2_pi_inv{0.918938533204673};
  auto log_dnrm = [&](double const x){
    return -x * x / 2. - log_sqrt_2_pi_inv;
  };

  bool const f_ub = std::isinf(ub),
             f_lb = std::isinf(lw);

  double const p_ub = f_ub ? 1 : pnorm_std(ub, 1L, 0L),
               p_lb = f_lb ? 0 : pnorm_std(lw, 1L, 0L),
               d_ub = f_ub ? 0 : std::exp(log_dnrm(ub) - pnorm_std(ub, 1L, 1L)),
               d_lb = f_lb ? 0 : std::exp(log_dnrm(lw) - pnorm_std(-lw, 1L, 1L)),
               d_ub_ub = f_ub ? 0 : ub * d_ub,
               d_lb_lb = f_lb ? 0 : lw * d_lb,
               sd_inv = 1 / *Sig_chol_tri();

  out[0] = p_ub - p_lb;
  out[1] = -(d_ub - d_lb) * sd_inv;
  out[2] = -(d_ub_ub - d_lb_lb) / 2 * sd_inv * sd_inv;
}

generic_l_factor::out_type generic_l_factor::get_output
  (double * res,  double const * sdest, size_t const minvls, int const inform,
   double const abserr, int const *indices){
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


cache_mem<double> likelihood::dmen;
cache_mem<double> pedigree_l_factor::dmem;
cache_mem<double> pedigree_l_factor_Hessian::dmem;
cache_mem<double> generic_l_factor::dmem;
cache_mem<int> pedigree_l_factor::imem;

template class cdf<pedigree_l_factor_Hessian>;

} // namespace pedmod
