// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mvndst
Rcpp::NumericVector mvndst(arma::vec const& lower, arma::vec const& upper, arma::vec const& mu, arma::mat const& sigma, unsigned const maxvls, double const abs_eps, double const rel_eps, int minvls, bool const do_reorder, bool const use_aprx);
RcppExport SEXP _pedmod_mvndst(SEXP lowerSEXP, SEXP upperSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP maxvlsSEXP, SEXP abs_epsSEXP, SEXP rel_epsSEXP, SEXP minvlsSEXP, SEXP do_reorderSEXP, SEXP use_aprxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec const& >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< unsigned const >::type maxvls(maxvlsSEXP);
    Rcpp::traits::input_parameter< double const >::type abs_eps(abs_epsSEXP);
    Rcpp::traits::input_parameter< double const >::type rel_eps(rel_epsSEXP);
    Rcpp::traits::input_parameter< int >::type minvls(minvlsSEXP);
    Rcpp::traits::input_parameter< bool const >::type do_reorder(do_reorderSEXP);
    Rcpp::traits::input_parameter< bool const >::type use_aprx(use_aprxSEXP);
    rcpp_result_gen = Rcpp::wrap(mvndst(lower, upper, mu, sigma, maxvls, abs_eps, rel_eps, minvls, do_reorder, use_aprx));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigree_ll_terms
SEXP get_pedigree_ll_terms(Rcpp::List data, unsigned const max_threads);
RcppExport SEXP _pedmod_get_pedigree_ll_terms(SEXP dataSEXP, SEXP max_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< unsigned const >::type max_threads(max_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree_ll_terms(data, max_threads));
    return rcpp_result_gen;
END_RCPP
}
// get_n_scales
int get_n_scales(SEXP ptr);
RcppExport SEXP _pedmod_get_n_scales(SEXP ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ptr(ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(get_n_scales(ptr));
    return rcpp_result_gen;
END_RCPP
}
// get_n_terms
int get_n_terms(SEXP ptr);
RcppExport SEXP _pedmod_get_n_terms(SEXP ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ptr(ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(get_n_terms(ptr));
    return rcpp_result_gen;
END_RCPP
}
// eval_pedigree_ll
Rcpp::NumericVector eval_pedigree_ll(SEXP ptr, arma::vec par, int const maxvls, double const abs_eps, double const rel_eps, Rcpp::Nullable<Rcpp::IntegerVector> indices, int const minvls, bool const do_reorder, bool const use_aprx, unsigned n_threads, Rcpp::Nullable<Rcpp::NumericVector> cluster_weights);
RcppExport SEXP _pedmod_eval_pedigree_ll(SEXP ptrSEXP, SEXP parSEXP, SEXP maxvlsSEXP, SEXP abs_epsSEXP, SEXP rel_epsSEXP, SEXP indicesSEXP, SEXP minvlsSEXP, SEXP do_reorderSEXP, SEXP use_aprxSEXP, SEXP n_threadsSEXP, SEXP cluster_weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ptr(ptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< int const >::type maxvls(maxvlsSEXP);
    Rcpp::traits::input_parameter< double const >::type abs_eps(abs_epsSEXP);
    Rcpp::traits::input_parameter< double const >::type rel_eps(rel_epsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerVector> >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< int const >::type minvls(minvlsSEXP);
    Rcpp::traits::input_parameter< bool const >::type do_reorder(do_reorderSEXP);
    Rcpp::traits::input_parameter< bool const >::type use_aprx(use_aprxSEXP);
    Rcpp::traits::input_parameter< unsigned >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type cluster_weights(cluster_weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(eval_pedigree_ll(ptr, par, maxvls, abs_eps, rel_eps, indices, minvls, do_reorder, use_aprx, n_threads, cluster_weights));
    return rcpp_result_gen;
END_RCPP
}
// eval_pedigree_grad
Rcpp::NumericVector eval_pedigree_grad(SEXP ptr, arma::vec par, int const maxvls, double const abs_eps, double const rel_eps, Rcpp::Nullable<Rcpp::IntegerVector> indices, int const minvls, bool const do_reorder, bool const use_aprx, unsigned n_threads, Rcpp::Nullable<Rcpp::NumericVector> cluster_weights);
RcppExport SEXP _pedmod_eval_pedigree_grad(SEXP ptrSEXP, SEXP parSEXP, SEXP maxvlsSEXP, SEXP abs_epsSEXP, SEXP rel_epsSEXP, SEXP indicesSEXP, SEXP minvlsSEXP, SEXP do_reorderSEXP, SEXP use_aprxSEXP, SEXP n_threadsSEXP, SEXP cluster_weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ptr(ptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< int const >::type maxvls(maxvlsSEXP);
    Rcpp::traits::input_parameter< double const >::type abs_eps(abs_epsSEXP);
    Rcpp::traits::input_parameter< double const >::type rel_eps(rel_epsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerVector> >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< int const >::type minvls(minvlsSEXP);
    Rcpp::traits::input_parameter< bool const >::type do_reorder(do_reorderSEXP);
    Rcpp::traits::input_parameter< bool const >::type use_aprx(use_aprxSEXP);
    Rcpp::traits::input_parameter< unsigned >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type cluster_weights(cluster_weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(eval_pedigree_grad(ptr, par, maxvls, abs_eps, rel_eps, indices, minvls, do_reorder, use_aprx, n_threads, cluster_weights));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP run_testthat_tests(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_pedmod_mvndst", (DL_FUNC) &_pedmod_mvndst, 10},
    {"_pedmod_get_pedigree_ll_terms", (DL_FUNC) &_pedmod_get_pedigree_ll_terms, 2},
    {"_pedmod_get_n_scales", (DL_FUNC) &_pedmod_get_n_scales, 1},
    {"_pedmod_get_n_terms", (DL_FUNC) &_pedmod_get_n_terms, 1},
    {"_pedmod_eval_pedigree_ll", (DL_FUNC) &_pedmod_eval_pedigree_ll, 11},
    {"_pedmod_eval_pedigree_grad", (DL_FUNC) &_pedmod_eval_pedigree_grad, 11},
    {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_pedmod(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
