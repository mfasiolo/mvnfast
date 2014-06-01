#include "mvnfast.h"
#include "internal.h"

/*
 *  Fast computation of Mahalanobis distance
 *
 * See ?maha() for a description of the arguments and output.
 */

/*
 * Interface to R
 */
SEXP mahaCpp(SEXP X, SEXP mu, SEXP sigma, SEXP ncores, SEXP isChol)
{
    using namespace Rcpp;
    
    try{
      
      arma::mat X_ = as<arma::mat>(X);
      arma::vec mu_ = as<arma::vec>(mu);  
      arma::mat sigma_ = as<arma::mat>(sigma); 
      unsigned int ncores_ = as<unsigned int>(ncores);
      bool isChol_ = as<bool>(isChol);
        
      NumericVector dist = wrap( mahaInt(X_, mu_, sigma_, ncores_, isChol_) );
      dist.attr( "dim" ) = R_NilValue;
      
      return dist;
      
    } catch( std::exception& __ex__){
      forward_exception_to_r(__ex__);
    } catch(...){
      ::Rf_error( "c++ exception (unknown reason)" );
    }
    
}
  

/* 
 *  Internal C++ function for Mahalanobis distance
*/
arma::vec mahaInt(arma::mat & X,  
                  arma::vec & mu,  
                  arma::mat & sigma,
                  const unsigned int ncores,
                  const bool isChol = false)
{
  using namespace arma;
  
  // Some sanity checks 
  if(mu.n_elem != sigma.n_cols) Rcpp::stop("The mean vector has a different dimensions from the covariance matrix.");
  if(X.n_cols != sigma.n_cols)  Rcpp::stop("The number of columns of X is different from the dimension of the covariance matrix.");
                   
  // Calculate transposed cholesky dec. unless sigma is alread a cholesky dec.
  mat cholDec;
  if( isChol == false ) {
     cholDec = trimatl(chol(sigma).t());
  }
  else{
     cholDec = trimatl(sigma.t()); 
     if(any(cholDec.diag() <= 0.0))  Rcpp::stop("The supplied cholesky decomposition has values <= 0.0 on the main diagonal.");
  }
  
  vec D = cholDec.diag();
    
  vec out(X.n_rows);
  
  #ifdef SUPPORT_OPENMP
  #pragma omp parallel num_threads(ncores) if(ncores > 1)                       
  {
  #endif
  
  // Declaring some private variables
  uint32_t d = X.n_cols;
  uint32_t n = X.n_rows;
  
  vec tmp(d);  
    
  double acc;
  uint32_t icol, irow, ii;  
  
  // For each of the "n" random vectors, forwardsolve the corresponding linear system.
  // Forwardsolve because I'm using the lower triangle Cholesky.
  #ifdef SUPPORT_OPENMP
  #pragma omp for schedule(static)
  #endif
  for(icol = 0; icol < n; icol++)
  {
        
    for(irow = 0; irow < d; irow++)
    {
     acc = 0.0;
     
     for(ii = 0; ii < irow; ii++) acc += tmp.at(ii) * cholDec.at(irow, ii);
     
     tmp.at(irow) = ( X.at(icol, irow) - mu.at(irow) - acc ) / D.at(irow);
    }
    
    out.at(icol) = sum(square(tmp)); 
  }
  
  #ifdef SUPPORT_OPENMP
  }
  #endif
    
  return out;
}




/* 
  #Equivalent R function:
  
  .fastMahalanobis <- function(X, mean, mcov)
  {
    dec <- chol(mcov)
    tmp <- forwardsolve(t(dec), t(X) - mean )
    colSums( tmp ^ 2 )
  }

*/
