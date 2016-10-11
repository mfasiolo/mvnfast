#ifndef _MVNFAST_H
#define _MVNFAST_H

#define BOOST_DISABLE_ASSERTS true

#include <RcppArmadillo.h>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/chi_squared_distribution.hpp>
#include "sitmo/sitmo.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

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
 
arma::vec mahaInt(arma::mat & X, arma::vec & mu, arma::mat & sigma, unsigned int ncores, bool isChol);

arma::vec dmvtInt( arma::mat X, arma::vec mu, arma::mat cholDec, bool log, double df, unsigned int ncores);

#endif
