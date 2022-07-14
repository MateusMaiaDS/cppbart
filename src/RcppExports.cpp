// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dhcauchy
double dhcauchy(double x, double location, double sigma);
RcppExport SEXP _cppbart_dhcauchy(SEXP xSEXP, SEXP locationSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type location(locationSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(dhcauchy(x, location, sigma));
    return rcpp_result_gen;
END_RCPP
}
// bart
List bart(Eigen::MatrixXd x_train, Eigen::VectorXd y, Eigen::MatrixXd x_test, int n_tree, int n_mcmc, int n_burn, int n_min_size, double tau, double mu, double tau_mu, double naive_sigma, double a_tau, double d_tau, double alpha, double beta);
RcppExport SEXP _cppbart_bart(SEXP x_trainSEXP, SEXP ySEXP, SEXP x_testSEXP, SEXP n_treeSEXP, SEXP n_mcmcSEXP, SEXP n_burnSEXP, SEXP n_min_sizeSEXP, SEXP tauSEXP, SEXP muSEXP, SEXP tau_muSEXP, SEXP naive_sigmaSEXP, SEXP a_tauSEXP, SEXP d_tauSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x_train(x_trainSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x_test(x_testSEXP);
    Rcpp::traits::input_parameter< int >::type n_tree(n_treeSEXP);
    Rcpp::traits::input_parameter< int >::type n_mcmc(n_mcmcSEXP);
    Rcpp::traits::input_parameter< int >::type n_burn(n_burnSEXP);
    Rcpp::traits::input_parameter< int >::type n_min_size(n_min_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type tau_mu(tau_muSEXP);
    Rcpp::traits::input_parameter< double >::type naive_sigma(naive_sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type a_tau(a_tauSEXP);
    Rcpp::traits::input_parameter< double >::type d_tau(d_tauSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(bart(x_train, y, x_test, n_tree, n_mcmc, n_burn, n_min_size, tau, mu, tau_mu, naive_sigma, a_tau, d_tau, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cppbart_dhcauchy", (DL_FUNC) &_cppbart_dhcauchy, 3},
    {"_cppbart_bart", (DL_FUNC) &_cppbart_bart, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_cppbart(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
