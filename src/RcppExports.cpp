// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// get_biconnected_components
Rcpp::List get_biconnected_components(Rcpp::IntegerVector const from, Rcpp::IntegerVector const to, Rcpp::IntegerVector const weights_ids, Rcpp::NumericVector const weights, Rcpp::NumericVector const edge_weights);
RcppExport SEXP _pedmod_get_biconnected_components(SEXP fromSEXP, SEXP toSEXP, SEXP weights_idsSEXP, SEXP weightsSEXP, SEXP edge_weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector const >::type from(fromSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector const >::type to(toSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector const >::type weights_ids(weights_idsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector const >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector const >::type edge_weights(edge_weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_biconnected_components(from, to, weights_ids, weights, edge_weights));
    return rcpp_result_gen;
END_RCPP
}
// get_block_cut_tree
Rcpp::List get_block_cut_tree(Rcpp::IntegerVector const from, Rcpp::IntegerVector const to, Rcpp::IntegerVector const weights_ids, Rcpp::NumericVector const weights, Rcpp::NumericVector const edge_weights);
RcppExport SEXP _pedmod_get_block_cut_tree(SEXP fromSEXP, SEXP toSEXP, SEXP weights_idsSEXP, SEXP weightsSEXP, SEXP edge_weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector const >::type from(fromSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector const >::type to(toSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector const >::type weights_ids(weights_idsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector const >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector const >::type edge_weights(edge_weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_block_cut_tree(from, to, weights_ids, weights, edge_weights));
    return rcpp_result_gen;
END_RCPP
}
// get_max_balanced_partition
Rcpp::List get_max_balanced_partition(Rcpp::IntegerVector const from, Rcpp::IntegerVector const to, Rcpp::IntegerVector const weights_ids, Rcpp::NumericVector const weights, Rcpp::NumericVector const edge_weights, double const slack, unsigned const max_kl_it_inner, unsigned const max_kl_it, unsigned const trace, bool const check_weights, bool const do_reorder);
RcppExport SEXP _pedmod_get_max_balanced_partition(SEXP fromSEXP, SEXP toSEXP, SEXP weights_idsSEXP, SEXP weightsSEXP, SEXP edge_weightsSEXP, SEXP slackSEXP, SEXP max_kl_it_innerSEXP, SEXP max_kl_itSEXP, SEXP traceSEXP, SEXP check_weightsSEXP, SEXP do_reorderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector const >::type from(fromSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector const >::type to(toSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector const >::type weights_ids(weights_idsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector const >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector const >::type edge_weights(edge_weightsSEXP);
    Rcpp::traits::input_parameter< double const >::type slack(slackSEXP);
    Rcpp::traits::input_parameter< unsigned const >::type max_kl_it_inner(max_kl_it_innerSEXP);
    Rcpp::traits::input_parameter< unsigned const >::type max_kl_it(max_kl_itSEXP);
    Rcpp::traits::input_parameter< unsigned const >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< bool const >::type check_weights(check_weightsSEXP);
    Rcpp::traits::input_parameter< bool const >::type do_reorder(do_reorderSEXP);
    rcpp_result_gen = Rcpp::wrap(get_max_balanced_partition(from, to, weights_ids, weights, edge_weights, slack, max_kl_it_inner, max_kl_it, trace, check_weights, do_reorder));
    return rcpp_result_gen;
END_RCPP
}
// unconnected_partition_rcpp
Rcpp::List unconnected_partition_rcpp(Rcpp::IntegerVector const from, Rcpp::IntegerVector const to, Rcpp::IntegerVector const weights_ids, Rcpp::NumericVector const weights, Rcpp::NumericVector const edge_weights, double const slack, unsigned const max_kl_it_inner, unsigned const max_kl_it, unsigned const trace, Rcpp::IntegerVector const init);
RcppExport SEXP _pedmod_unconnected_partition_rcpp(SEXP fromSEXP, SEXP toSEXP, SEXP weights_idsSEXP, SEXP weightsSEXP, SEXP edge_weightsSEXP, SEXP slackSEXP, SEXP max_kl_it_innerSEXP, SEXP max_kl_itSEXP, SEXP traceSEXP, SEXP initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector const >::type from(fromSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector const >::type to(toSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector const >::type weights_ids(weights_idsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector const >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector const >::type edge_weights(edge_weightsSEXP);
    Rcpp::traits::input_parameter< double const >::type slack(slackSEXP);
    Rcpp::traits::input_parameter< unsigned const >::type max_kl_it_inner(max_kl_it_innerSEXP);
    Rcpp::traits::input_parameter< unsigned const >::type max_kl_it(max_kl_itSEXP);
    Rcpp::traits::input_parameter< unsigned const >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector const >::type init(initSEXP);
    rcpp_result_gen = Rcpp::wrap(unconnected_partition_rcpp(from, to, weights_ids, weights, edge_weights, slack, max_kl_it_inner, max_kl_it, trace, init));
    return rcpp_result_gen;
END_RCPP
}
// mvndst
Rcpp::NumericVector mvndst(arma::vec const& lower, arma::vec const& upper, arma::vec const& mu, arma::mat const& sigma, unsigned const maxvls, double const abs_eps, double const rel_eps, int minvls, bool const do_reorder, bool const use_aprx, int const method, unsigned const n_sequences);
RcppExport SEXP _pedmod_mvndst(SEXP lowerSEXP, SEXP upperSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP maxvlsSEXP, SEXP abs_epsSEXP, SEXP rel_epsSEXP, SEXP minvlsSEXP, SEXP do_reorderSEXP, SEXP use_aprxSEXP, SEXP methodSEXP, SEXP n_sequencesSEXP) {
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
    Rcpp::traits::input_parameter< int const >::type method(methodSEXP);
    Rcpp::traits::input_parameter< unsigned const >::type n_sequences(n_sequencesSEXP);
    rcpp_result_gen = Rcpp::wrap(mvndst(lower, upper, mu, sigma, maxvls, abs_eps, rel_eps, minvls, do_reorder, use_aprx, method, n_sequences));
    return rcpp_result_gen;
END_RCPP
}
// get_pedigree_ll_terms
SEXP get_pedigree_ll_terms(Rcpp::List data, unsigned const max_threads, unsigned const n_sequences);
RcppExport SEXP _pedmod_get_pedigree_ll_terms(SEXP dataSEXP, SEXP max_threadsSEXP, SEXP n_sequencesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< unsigned const >::type max_threads(max_threadsSEXP);
    Rcpp::traits::input_parameter< unsigned const >::type n_sequences(n_sequencesSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pedigree_ll_terms(data, max_threads, n_sequences));
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
Rcpp::NumericVector eval_pedigree_ll(SEXP ptr, arma::vec par, int const maxvls, double const abs_eps, double const rel_eps, Rcpp::Nullable<Rcpp::IntegerVector> indices, int const minvls, bool const do_reorder, bool const use_aprx, unsigned n_threads, Rcpp::Nullable<Rcpp::NumericVector> cluster_weights, int const method);
RcppExport SEXP _pedmod_eval_pedigree_ll(SEXP ptrSEXP, SEXP parSEXP, SEXP maxvlsSEXP, SEXP abs_epsSEXP, SEXP rel_epsSEXP, SEXP indicesSEXP, SEXP minvlsSEXP, SEXP do_reorderSEXP, SEXP use_aprxSEXP, SEXP n_threadsSEXP, SEXP cluster_weightsSEXP, SEXP methodSEXP) {
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
    Rcpp::traits::input_parameter< int const >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(eval_pedigree_ll(ptr, par, maxvls, abs_eps, rel_eps, indices, minvls, do_reorder, use_aprx, n_threads, cluster_weights, method));
    return rcpp_result_gen;
END_RCPP
}
// eval_pedigree_grad
Rcpp::NumericVector eval_pedigree_grad(SEXP ptr, arma::vec par, int const maxvls, double const abs_eps, double const rel_eps, Rcpp::Nullable<Rcpp::IntegerVector> indices, int const minvls, bool const do_reorder, bool const use_aprx, unsigned n_threads, Rcpp::Nullable<Rcpp::NumericVector> cluster_weights, int const method);
RcppExport SEXP _pedmod_eval_pedigree_grad(SEXP ptrSEXP, SEXP parSEXP, SEXP maxvlsSEXP, SEXP abs_epsSEXP, SEXP rel_epsSEXP, SEXP indicesSEXP, SEXP minvlsSEXP, SEXP do_reorderSEXP, SEXP use_aprxSEXP, SEXP n_threadsSEXP, SEXP cluster_weightsSEXP, SEXP methodSEXP) {
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
    Rcpp::traits::input_parameter< int const >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(eval_pedigree_grad(ptr, par, maxvls, abs_eps, rel_eps, indices, minvls, do_reorder, use_aprx, n_threads, cluster_weights, method));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP run_testthat_tests(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_pedmod_get_biconnected_components", (DL_FUNC) &_pedmod_get_biconnected_components, 5},
    {"_pedmod_get_block_cut_tree", (DL_FUNC) &_pedmod_get_block_cut_tree, 5},
    {"_pedmod_get_max_balanced_partition", (DL_FUNC) &_pedmod_get_max_balanced_partition, 11},
    {"_pedmod_unconnected_partition_rcpp", (DL_FUNC) &_pedmod_unconnected_partition_rcpp, 10},
    {"_pedmod_mvndst", (DL_FUNC) &_pedmod_mvndst, 12},
    {"_pedmod_get_pedigree_ll_terms", (DL_FUNC) &_pedmod_get_pedigree_ll_terms, 3},
    {"_pedmod_get_n_scales", (DL_FUNC) &_pedmod_get_n_scales, 1},
    {"_pedmod_get_n_terms", (DL_FUNC) &_pedmod_get_n_terms, 1},
    {"_pedmod_eval_pedigree_ll", (DL_FUNC) &_pedmod_eval_pedigree_ll, 12},
    {"_pedmod_eval_pedigree_grad", (DL_FUNC) &_pedmod_eval_pedigree_grad, 12},
    {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_pedmod(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
