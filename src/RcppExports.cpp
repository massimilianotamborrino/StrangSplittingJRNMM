// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// mv_mult_JRNMM_
NumericVector mv_mult_JRNMM_(NumericMatrix mat, NumericVector vec);
RcppExport SEXP _StrangSplittingJRNMM_mv_mult_JRNMM_(SEXP matSEXP, SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(mv_mult_JRNMM_(mat, vec));
    return rcpp_result_gen;
END_RCPP
}
// sigmoid_JRNMM_Cpp_
double sigmoid_JRNMM_Cpp_(double x, double vmax, double v0, double r);
RcppExport SEXP _StrangSplittingJRNMM_sigmoid_JRNMM_Cpp_(SEXP xSEXP, SEXP vmaxSEXP, SEXP v0SEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type vmax(vmaxSEXP);
    Rcpp::traits::input_parameter< double >::type v0(v0SEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(sigmoid_JRNMM_Cpp_(x, vmax, v0, r));
    return rcpp_result_gen;
END_RCPP
}
// exponential_matrix_JRNMM_
NumericMatrix exponential_matrix_JRNMM_(NumericVector dGamma, double h);
RcppExport SEXP _StrangSplittingJRNMM_exponential_matrix_JRNMM_(SEXP dGammaSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type dGamma(dGammaSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(exponential_matrix_JRNMM_(dGamma, h));
    return rcpp_result_gen;
END_RCPP
}
// covariance_matrix_JRNMM_
NumericMatrix covariance_matrix_JRNMM_(NumericVector dGamma, NumericVector dSigma, double h);
RcppExport SEXP _StrangSplittingJRNMM_covariance_matrix_JRNMM_(SEXP dGammaSEXP, SEXP dSigmaSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type dGamma(dGammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dSigma(dSigmaSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(covariance_matrix_JRNMM_(dGamma, dSigma, h));
    return rcpp_result_gen;
END_RCPP
}
// linear_JRNMM_Cpp_
NumericVector linear_JRNMM_Cpp_(NumericVector vec, NumericMatrix dm, NumericVector xi);
RcppExport SEXP _StrangSplittingJRNMM_linear_JRNMM_Cpp_(SEXP vecSEXP, SEXP dmSEXP, SEXP xiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dm(dmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xi(xiSEXP);
    rcpp_result_gen = Rcpp::wrap(linear_JRNMM_Cpp_(vec, dm, xi));
    return rcpp_result_gen;
END_RCPP
}
// nonlinear_JRNMM_Cpp_
NumericVector nonlinear_JRNMM_Cpp_(int N, NumericVector vec, double h, NumericMatrix Theta, NumericMatrix Rho, NumericMatrix K);
RcppExport SEXP _StrangSplittingJRNMM_nonlinear_JRNMM_Cpp_(SEXP NSEXP, SEXP vecSEXP, SEXP hSEXP, SEXP ThetaSEXP, SEXP RhoSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Rho(RhoSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(nonlinear_JRNMM_Cpp_(N, vec, h, Theta, Rho, K));
    return rcpp_result_gen;
END_RCPP
}
// JRNMM_Splitting_Cpp_
NumericMatrix JRNMM_Splitting_Cpp_(int N_i, NumericVector grid_i, double h_i, NumericVector startv, NumericMatrix dm_i, NumericVector meanVec_i, NumericMatrix covMat_i, NumericMatrix Theta_i, NumericMatrix Rho_i, NumericMatrix K_i);
RcppExport SEXP _StrangSplittingJRNMM_JRNMM_Splitting_Cpp_(SEXP N_iSEXP, SEXP grid_iSEXP, SEXP h_iSEXP, SEXP startvSEXP, SEXP dm_iSEXP, SEXP meanVec_iSEXP, SEXP covMat_iSEXP, SEXP Theta_iSEXP, SEXP Rho_iSEXP, SEXP K_iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N_i(N_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type grid_i(grid_iSEXP);
    Rcpp::traits::input_parameter< double >::type h_i(h_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type startv(startvSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dm_i(dm_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type meanVec_i(meanVec_iSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covMat_i(covMat_iSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Theta_i(Theta_iSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Rho_i(Rho_iSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type K_i(K_iSEXP);
    rcpp_result_gen = Rcpp::wrap(JRNMM_Splitting_Cpp_(N_i, grid_i, h_i, startv, dm_i, meanVec_i, covMat_i, Theta_i, Rho_i, K_i));
    return rcpp_result_gen;
END_RCPP
}
// meanHO_
NumericVector meanHO_(int N, double h, NumericMatrix expM, NumericVector x);
RcppExport SEXP _StrangSplittingJRNMM_meanHO_(SEXP NSEXP, SEXP hSEXP, SEXP expMSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type expM(expMSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(meanHO_(N, h, expM, x));
    return rcpp_result_gen;
END_RCPP
}
// fastJRNMM_Splitting_Cpp_
NumericMatrix fastJRNMM_Splitting_Cpp_(int N, NumericVector grid, double h, NumericVector start, NumericVector dGamma, NumericVector dSigma, NumericMatrix Theta, NumericMatrix Rho, NumericMatrix K);
RcppExport SEXP _StrangSplittingJRNMM_fastJRNMM_Splitting_Cpp_(SEXP NSEXP, SEXP gridSEXP, SEXP hSEXP, SEXP startSEXP, SEXP dGammaSEXP, SEXP dSigmaSEXP, SEXP ThetaSEXP, SEXP RhoSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dGamma(dGammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dSigma(dSigmaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Rho(RhoSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(fastJRNMM_Splitting_Cpp_(N, grid, h, start, dGamma, dSigma, Theta, Rho, K));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_StrangSplittingJRNMM_mv_mult_JRNMM_", (DL_FUNC) &_StrangSplittingJRNMM_mv_mult_JRNMM_, 2},
    {"_StrangSplittingJRNMM_sigmoid_JRNMM_Cpp_", (DL_FUNC) &_StrangSplittingJRNMM_sigmoid_JRNMM_Cpp_, 4},
    {"_StrangSplittingJRNMM_exponential_matrix_JRNMM_", (DL_FUNC) &_StrangSplittingJRNMM_exponential_matrix_JRNMM_, 2},
    {"_StrangSplittingJRNMM_covariance_matrix_JRNMM_", (DL_FUNC) &_StrangSplittingJRNMM_covariance_matrix_JRNMM_, 3},
    {"_StrangSplittingJRNMM_linear_JRNMM_Cpp_", (DL_FUNC) &_StrangSplittingJRNMM_linear_JRNMM_Cpp_, 3},
    {"_StrangSplittingJRNMM_nonlinear_JRNMM_Cpp_", (DL_FUNC) &_StrangSplittingJRNMM_nonlinear_JRNMM_Cpp_, 6},
    {"_StrangSplittingJRNMM_JRNMM_Splitting_Cpp_", (DL_FUNC) &_StrangSplittingJRNMM_JRNMM_Splitting_Cpp_, 10},
    {"_StrangSplittingJRNMM_meanHO_", (DL_FUNC) &_StrangSplittingJRNMM_meanHO_, 4},
    {"_StrangSplittingJRNMM_fastJRNMM_Splitting_Cpp_", (DL_FUNC) &_StrangSplittingJRNMM_fastJRNMM_Splitting_Cpp_, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_StrangSplittingJRNMM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
