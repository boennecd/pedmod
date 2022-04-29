#include "cdfaprx.h"

namespace pedmod {

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
               d_lb = f_lb ? 0 : std::exp(log_dnrm(lw) - pnorm_std(lw, 1L, 1L)),
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
cache_mem<double> generic_l_factor::dmem;
cache_mem<int> pedigree_l_factor::imem;

} // namespace pedmod {
