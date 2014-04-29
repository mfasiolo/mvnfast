#include "mvnfast.h"
#include "internal.h"

/*
  Fast computation of pdf of a multivariate normal distribution
*/
  
  SEXP dmvnCpp(SEXP X_,  
               SEXP mu_,  
               SEXP sigma_, 
               SEXP log_,
               SEXP ncores_,
               SEXP isChol_) 
{ 
    using namespace arma;
    
    try{
      mat X = Rcpp::as<mat>(X_);
      vec mu = Rcpp::as<vec>(mu_);  
      mat sigma = Rcpp::as<mat>(sigma_); 
      bool log = Rcpp::as<bool>(log_); 
      int  ncores = Rcpp::as<int>(ncores_); 
      bool isChol = Rcpp::as<bool>(isChol_); 
      
      int d = X.n_cols;
      
      // Calculate cholesky dec. unless sigma is alread a cholesky dec.
      mat cholDec;
      if( isChol == false ) {
        cholDec = chol(sigma);
      }
      else{
        cholDec = sigma;
      }
      
      // Calculate residual sum of squares
      vec out = - 0.5 * mahaInt(X, mu, cholDec, ncores, true);
        
      out = out - ( (d / 2.0) * std::log(2.0 * M_PI) + sum(arma::log(cholDec.diag())) );
      
      if (log == false) {
        out = exp(out);
      }
      
      // Dropping the dimensionality of the output vector
      Rcpp::NumericVector Rout = Rcpp::wrap(out);
      Rout.attr( "dim" ) = R_NilValue;
      
      return Rout;
      
    } catch( std::exception& __ex__){
      forward_exception_to_r(__ex__);
    } catch(...){
      ::Rf_error( "c++ exception (unknown reason)" );
    }
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
