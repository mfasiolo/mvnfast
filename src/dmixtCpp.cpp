/* 
 Copyright (C) 2014 Matteo Fasiolo  matteo.fasiolo@gmail.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
(www.gnu.org/copyleft/gpl.html)

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
USA. */

#include "mvnfast.h"

/*
* Fast computation of pdf of a multivariate normal distribution
*
* See ?dmixt() for a description of the arguments and output.
*/

/*
* Interface to R
*/
RcppExport SEXP dmixtCpp(SEXP X_,  
                         SEXP mu_,  
                         SEXP sigma_, 
                         SEXP df_,
                         SEXP w_,
                         SEXP log_,
                         SEXP ncores_,
                         SEXP isChol_, 
                         SEXP A_) 
{ 
  using namespace Rcpp;
  
  try{
    arma::mat X = as<arma::mat>(X_);
    arma::mat mu = as<arma::mat>(mu_);  
    List sigma = List(sigma_);
    arma::vec w = as<arma::vec>(w_);
    double df = as<double>(df_); 
    bool log = as<bool>(log_); 
    unsigned int  ncores = as<unsigned int>(ncores_); 
    bool isChol = as<bool>(isChol_);
    NumericMatrix A = NumericMatrix(A_);
    
    // We A wrap into a arma::mat "tmp" without making a copy.
    arma::mat armA( A.begin(), A.nrow(), A.ncol(), false );
    
    unsigned int d = mu.n_cols; // Dimension of random vectors
    uint32_t m = mu.n_rows; // Number of mixture components
    uint32_t n = X.n_rows; // Number of samples
    
    if(ncores == 0) stop("ncores has to be positive");
    if(mu.n_rows != m) stop("mu.n_rows != m");
    if(w.n_elem != m) stop("w.n_elem != m");
    if( X.n_cols != d ) Rcpp::stop("X.n_cols != d"); 
    if( A.nrow() != n ) Rcpp::stop("A.nrow() != n");
    if( A.ncol() != m ) Rcpp::stop("A.ncol() != m");
    
    // Get list of Cholesky decompositions of covariance matrices
    arma::mat tmpMat;
    std::vector< arma::mat > cholDec;
    for(int ii = 0; ii < m; ii++)
    {
      tmpMat = as<arma::mat>(wrap(sigma[ii]));
      
      if(d != tmpMat.n_cols) stop("mu.n_cols != sigma[ii].n_cols");
      if(d != tmpMat.n_rows) stop("mu.n_cols != sigma[ii].n_rows");
      
      // Calculate cholesky dec unless sigma is already a cholesky dec.
      if( isChol == false ) {
        cholDec.push_back( trimatu(arma::chol(tmpMat)) );
      }
      else{
        cholDec.push_back( trimatu(tmpMat) );
      }
    }
    
    // Evaluate log-density of each mixture density over all samples.
    for(int imix = 0; imix<m; imix++)
    {
      armA.col(imix) = dmvtInt(X, mu.row(imix).t(), cholDec[imix], true, df, ncores);
    }
    
    // Calculating density of each sample using log sum exp trick
    arma::vec out(n);
    double mx;
    for(int ii = 0; ii<n; ii++)
    {
      mx = armA.row(ii).max();
      armA.row(ii) -= mx;
      out.at(ii) = mx; 
    }
    
    // max(log(xi)) + log( sum(wi * exp(log(xi) - max(log(xi)))) )
    out = out + arma::log(exp(armA) * w);
    
    if( log == false ){ out = exp(out); }
    
    // Dropping the dimensionality of the output vector
    Rcpp::NumericVector Rout = Rcpp::wrap( out );
    Rout.attr( "dim" ) = R_NilValue;
    
    return Rout;
    
  } catch( std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return wrap(NA_REAL);
}

