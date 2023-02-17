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
  * See ?dmvt() for a description of the arguments and output.
*/
  
  /*
  * Interface to R
*/
  RcppExport SEXP dmvtCpp(SEXP X_,  
                          SEXP mu_,  
                          SEXP sigma_, 
                          SEXP df_, 
                          SEXP log_,
                          SEXP ncores_,
                          SEXP isChol_) 
{ 
  using namespace Rcpp;
  
  try{
    arma::mat X = as<arma::mat>(X_);
    arma::vec mu = as<arma::vec>(mu_);  
    arma::mat sigma = as<arma::mat>(sigma_); 
    double df = as<double>(df_); 
    bool log = as<bool>(log_); 
    unsigned int  ncores = as<unsigned int>(ncores_); 
    bool isChol = as<bool>(isChol_);
    
    if(ncores == 0) stop("ncores has to be positive.");
    if( X.n_cols != mu.n_elem ) Rcpp::stop("X.n_cols != mu.n_elem"); 
    if( X.n_cols != sigma.n_cols ) Rcpp::stop("X.n_cols != sigma.n_cols"); 
    if( sigma.n_rows != sigma.n_cols ) Rcpp::stop("sigma.n_rows != sigma.n_cols"); 
    
    #ifdef _OPENMP
    omp_set_num_threads(ncores);
    #endif
    
    // Calculate cholesky dec. unless sigma is alread a cholesky dec.
    arma::mat cholDec;
    if( isChol == false ) {
      cholDec = arma::chol(sigma);
    }
    else{
      cholDec = sigma;
    }
    
    // Dropping the dimensionality of the output vector
    Rcpp::NumericVector Rout = Rcpp::wrap( dmvtInt( X, mu, cholDec, log, df, ncores) );
    Rout.attr( "dim" ) = R_NilValue;
    
    return Rout;
    
  } catch( std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return wrap(NA_REAL);
  }


/* 
  *  Internal C++ function
*/
  arma::vec dmvtInt( arma::mat X, arma::vec mu, arma::mat cholDec, bool log, double df, unsigned int ncores)
{
  using namespace arma;
  
  unsigned int d = X.n_cols;
  
  vec out = mahaInt(X, mu, cholDec, ncores, true);
  
  if( df <= 0.0 ){ // Multivariate normal density OR ...
    
    out = - 0.5 * out - ( (d / 2.0) * std::log(2.0 * M_PI) + sum(arma::log(cholDec.diag())) );
    
  } else { // ... multivariate Student-t density
    
  #ifdef _OPENMP
  #pragma omp parallel num_threads(ncores) if(ncores > 1)
  {
  #endif
  
  uint32_t ii;  
  uint32_t n = X.n_rows;  
  double logDet = sum(arma::log(cholDec.diag())); 
  double c = lgamma((d + df)/2.0) - (lgamma(df/2.0) + logDet + d/2.0 * std::log(M_PI * df));

  #ifdef _OPENMP
  #pragma omp for schedule(static)
  #endif
  for(ii = 0; ii < n; ii++)
  {
     out.at(ii) = c - 0.5 * (df + d) * log1p(out.at(ii)/df);
  }
    
  #ifdef _OPENMP
  }
  #endif
    
  }
  
  if (log == false) out = exp(out);
  
  return( out );
}


/* 
  # Equivalent R function:
  
  myNorm <- function(X, mean, varcov)
  {
    dec <- chol(varcov)
    tmp <- forwardsolve(dec, t(X) - mean )  
    rss <- colSums( tmp ^ 2 )
    
    return( exp( - sum(log(diag(dec))) - 0.5 * length(mean) * log(2 * pi) - 0.5 * rss )  )
  }
*/
  