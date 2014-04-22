#include "mvn.h"
#include "internal.h"

/*
  * Fast computation of Mahalanobis distance
*/
  
  /*
  *  Interface with R
*/
  
  SEXP mahaCpp(SEXP X, SEXP mu, SEXP sigma, SEXP isChol)
{
    using namespace Rcpp;
    
    try{
      
      arma::mat X_ = as<arma::mat>(X);
      arma::vec mu_ = as<arma::vec>(mu);  
      arma::mat sigma_ = as<arma::mat>(sigma); 
      bool isChol_ = as<bool>(isChol);
      
      arma::rowvec dist = mahaInt(X_, mu_, sigma_, isChol_);
      
      return Rcpp::wrap(dist);
      
    } catch( std::exception& __ex__){
      forward_exception_to_r(__ex__);
    } catch(...){
      ::Rf_error( "c++ exception (unknown reason)" );
    }
    
  }


/* 
  *  Internal C++ function
*/
  
  arma::rowvec mahaInt(arma::mat & X,  
                       arma::vec & mu,  
                       arma::mat & sigma, 
                       bool isChol = false)
{
    using namespace arma;
    
    // n vectors of dimension d.
    int n = X.n_rows;
    int d = X.n_cols;
    
    // subtract mean from each row.
    X.each_row() -= mu.t();
    
    // Calculate transposed cholesky dec. unless sigma is alread a cholesky dec.
    mat cholDec;
    if( isChol == false ) {
      cholDec = trans( chol(sigma) );
    }
    else{
      cholDec = trans(sigma);
    }
    
    // Solving n linear system, result are normalized residuals.
    mat res = solve( trimatl( cholDec ), X.t() );
    
    // Residual sums of squares are the squared norm of each column of res. 
    rowvec rss = sum(square(res), 0);
    
    return rss;
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
