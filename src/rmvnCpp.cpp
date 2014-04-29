#include "mvnfast.h"
#include "internal.h"

/*
  Fast computation of pdf of a multivariate normal distribution
*/

SEXP rmvnCpp(SEXP n_,  
             SEXP mu_,  
             SEXP sigma_,
             SEXP ncores_,
             SEXP isChol_) 
{ 
    using namespace Rcpp;
    
    try{
      
      RNGScope scope;
      
      int n = as<int>(n_);
      arma::rowvec mu = as<arma::rowvec>(mu_);  
      arma::mat sigma = as<arma::mat>(sigma_); 
      int  ncores = as<int>(ncores_); 
      bool isChol = as<bool>(isChol_); 
      
      int d = mu.n_elem;
    
      // Calculate cholesky dec. unless sigma is alread a cholesky dec.
      arma::mat cholDec;
      if( isChol == false ) {
        cholDec = arma::chol(sigma);
      }
      else{
        cholDec = sigma;
      }
                  
      // Generate standard normals using R rng, put it into the tmp matrix, which is then
      // converted to a arma::mat, without copying.
      NumericMatrix out(n, d);
      arma::mat tmp( out.begin(), out.nrow(), out.ncol(), false );
      
      NumericVector seeds = runif(4, 1.0, 1.844674e+19);
      
      #pragma omp parallel num_threads(ncores) if(ncores > 1)
      {
      
      double acc;
      int irow, icol, ii;
      arma::rowvec work(d);
      
      std::mt19937_64 engine( static_cast<uint64_t>(seeds[omp_get_thread_num()]) );
      std::normal_distribution<> normal(0.0, 1.0);
      
      #pragma omp for schedule(static)
      for (int irow = 0; irow < n; irow++) 
        for (int icol = 0; icol < d; icol++) 
           out(irow, icol) = normal(engine);
        
      #pragma omp for schedule(static)
      for(irow = 0; irow < n; irow++)
      {
       
       for(icol = 0; icol < d; icol++)
       {
        acc = 0.0;
        
        for(ii = 0; ii <= icol; ii++) acc += tmp.at(irow, ii) * cholDec.at(ii, icol);
        
        work.at(icol) = acc;
        
       }
       
       tmp(arma::span(irow), arma::span::all) = work + mu;       
      }
      }
      
      return out;
            
    } catch( std::exception& __ex__){
      forward_exception_to_r(__ex__);
    } catch(...){
      ::Rf_error( "c++ exception (unknown reason)" );
    }
}