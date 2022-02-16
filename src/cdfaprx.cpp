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

cache_mem<double> likelihood::dmen;
cache_mem<double> pedigree_l_factor::dmem;
cache_mem<double> generic_l_factor::dmem;
cache_mem<int> pedigree_l_factor::imem;

} // namespace pedmod {
