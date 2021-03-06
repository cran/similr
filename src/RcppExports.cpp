// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// findThem
Rcpp::List findThem(std::vector<std::string> sources, std::vector<std::string> targets, int tophits, double min_similarity);
RcppExport SEXP _similr_findThem(SEXP sourcesSEXP, SEXP targetsSEXP, SEXP tophitsSEXP, SEXP min_similaritySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type sources(sourcesSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type targets(targetsSEXP);
    Rcpp::traits::input_parameter< int >::type tophits(tophitsSEXP);
    Rcpp::traits::input_parameter< double >::type min_similarity(min_similaritySEXP);
    rcpp_result_gen = Rcpp::wrap(findThem(sources, targets, tophits, min_similarity));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_similr_findThem", (DL_FUNC) &_similr_findThem, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_similr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
