#ifndef _SYNLIK_H
#define _SYNLIK_H

#include <RcppArmadillo.h>

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that 
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the 
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */
 
RcppExport SEXP cleanStats(SEXP inMat);

RcppExport SEXP checkBoundsCpp(SEXP theMean_, SEXP cholFact_, SEXP toCheck_, SEXP upper_, SEXP lower_, SEXP output_);

RcppExport SEXP mahaCpp(SEXP X, SEXP mu, SEXP sigma, SEXP isChol);

RcppExport SEXP dmvnCpp(SEXP X_, SEXP mu_, SEXP sigma_, SEXP log_, SEXP isChol_);

/*
 * Simulators
 */

RcppExport SEXP simpleModelsWrap(SEXP model, SEXP days, SEXP nSimul, SEXP params, SEXP nBurn, SEXP randInit, SEXP initVal);

#endif
