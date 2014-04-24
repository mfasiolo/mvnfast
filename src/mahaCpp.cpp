#include "mvn.h"
#include "internal.h"
#include <omp.h>

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
 *  Forward solve linear system
*/
arma::mat forwardSolve(arma::mat & A, arma::mat & X)
{
  using namespace arma;
  
  //omp_set_num_threads(2);
  
  int d = X.n_rows;
  int n = X.n_cols;
  mat out(d, n);
  
  //#pragma omp parallel for schedule(static)
  for(int ii = 0; ii < n; ii++)
  {
    out(span::all, ii) = solve(trimatu(A), X(span::all, ii));
  }
  
  return(out); 
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
      cholDec = chol(sigma);
    }
    else{
      cholDec = trans(sigma);
    }
    
    // Solving n linear system, result are normalized residuals.
    //X = X.t();
    //mat res = forwardSolve( cholDec, X );
    mat res = solve( trimatu( cholDec ), X.t() );
    
    // Residual sums of squares are the squared norm of each column of res. 
    rowvec rss = sum(square(res), 0.0);
    
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
