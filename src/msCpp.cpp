#include "mvnfast.h"
#include "internal.h"

/*
 *  Mean-shift algorithm
 *
 * See ?ms() for a description of the arguments and output.
*/

SEXP msCpp(SEXP init_, SEXP X_, SEXP cholDec_, SEXP ncores_, SEXP tol_, SEXP store_)
{
    using namespace Rcpp;
    
    try{
      
      arma::vec init = as<arma::vec>(init_);
      arma::mat X = as<arma::mat>(X_); 
      arma::mat cholDec = as<arma::mat>(cholDec_); 
      int ncores = as<int>(ncores_);
      double tol = as<double>(tol_);
      bool store = as<bool>(store_);
      
      int d = init.n_elem;
      int n = X.n_rows;
      
      if( d != X.n_cols ) stop( "The ncol(X) has to equal to length(init)" );
      
      arma::vec currPos = init; 
      arma::vec oldPos(d);
      arma::vec weights(n);
      std::list< std::vector<double> > traj;
            
      if(store) traj.push_back( as<std::vector<double>>( wrap(init) ) );
      
      do{
        
        oldPos = currPos; 
        weights = dmvnInt(X, oldPos, cholDec, false, ncores);
        currPos = arma::conv_to<arma::vec>::from( weights.t() * X / sum(weights) ); 
    
        if(store) traj.push_back( as<std::vector<double>>( wrap(currPos) ) );
        
      }
      while( any( abs(currPos - oldPos) > tol ) );
            
      return List::create( _["final"] = wrap( currPos ), _["traj"] = wrap( traj ) );
      
    } catch( std::exception& __ex__){
      forward_exception_to_r(__ex__);
    } catch(...){
      ::Rf_error( "c++ exception (unknown reason)" );
    }
    
  }