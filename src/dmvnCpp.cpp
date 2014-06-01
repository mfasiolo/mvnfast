#include "mvnfast.h"

/*
 * Fast computation of pdf of a multivariate normal distribution
 *
 * See ?dmvn() for a description of the arguments and output.
 */
 
/*
 * Interface to R
 */
RcppExport SEXP dmvnCpp(SEXP X_,  
             SEXP mu_,  
             SEXP sigma_, 
             SEXP log_,
             SEXP ncores_,
             SEXP isChol_) 
{ 
    using namespace Rcpp;
    
    try{
      arma::mat X = as<arma::mat>(X_);
      arma::vec mu = as<arma::vec>(mu_);  
      arma::mat sigma = as<arma::mat>(sigma_); 
      bool log = as<bool>(log_); 
      unsigned int  ncores = as<unsigned int>(ncores_); 
      bool isChol = as<bool>(isChol_);
      
      // Calculate cholesky dec. unless sigma is alread a cholesky dec.
     arma::mat cholDec;
     if( isChol == false ) {
        cholDec = arma::chol(sigma);
     }
     else{
        cholDec = sigma;
     }
      
      // Dropping the dimensionality of the output vector
      Rcpp::NumericVector Rout = Rcpp::wrap( dmvnInt( X, mu, cholDec, log, ncores) );
      Rout.attr( "dim" ) = R_NilValue;
      
      return Rout;
      
    } catch( std::exception& __ex__){
      forward_exception_to_r(__ex__);
    } catch(...){
      ::Rf_error( "c++ exception (unknown reason)" );
    }
}


/* 
 *  Internal C++ function
*/
arma::vec dmvnInt( arma::mat X, arma::vec mu, arma::mat cholDec, bool log, unsigned int ncores)
{
 using namespace arma;
  
 unsigned int d = X.n_cols;
      
 vec out = - 0.5 * mahaInt(X, mu, cholDec, ncores, true);
        
 out = out - ( (d / 2.0) * std::log(2.0 * M_PI) + sum(arma::log(cholDec.diag())) );
      
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
