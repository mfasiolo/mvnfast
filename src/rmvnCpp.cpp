#include "mvnfast.h"

/*
 * Simulate random variables from a multivariate normal distribution
 *
 * See ?rmvn() for a description of the arguments and output.
*/

RcppExport SEXP rmvnCpp(SEXP n_,  
             SEXP mu_,  
             SEXP sigma_,
             SEXP ncores_,
             SEXP isChol_) 
{ 
    using namespace Rcpp;
    
    try{
      
      RNGScope scope;
      
      uint32_t n = as<uint32_t>(n_);
      arma::rowvec mu = as<arma::rowvec>(mu_);  
      arma::mat sigma = as<arma::mat>(sigma_); 
      unsigned int  ncores = as<unsigned int>(ncores_); 
      bool isChol = as<bool>(isChol_); 
      
      unsigned int d = mu.n_elem;
    
      // Calculate cholesky dec unless sigma is already a cholesky dec.
      arma::mat cholDec;
      if( isChol == false ) {
        cholDec = trimatu( arma::chol(sigma) );
      }
      else{
        cholDec = trimatu( sigma );
      }
                  
      // This "out" the matrix that will be filled with firstly with standard normal rvs,
      // and finally with multivariate normal rvs.
      // We wrap into a arma::mat "tmp" without making a copy.
      NumericMatrix out(n, d);
      arma::mat tmp( out.begin(), out.nrow(), out.ncol(), false );
      
      // What I do to seed the sitmo::prng_engine is tricky. I produce "ncores" uniform numbers between 1 and the largest uint32_t,
      // which are the seeds. I put the first one in "coreSeed". If there is no support for OpenMP only this seed
      // will be used, as the computations will be sequential. If there is support for OpenMP, "coreSeed" will
      // be over-written, so that each core will get its own seed.
      NumericVector seeds = runif(ncores, 1.0, std::numeric_limits<uint32_t>::max());
      
      #ifdef SUPPORT_OPENMP
      #pragma omp parallel num_threads(ncores) if(ncores > 1)
      {
      #endif
      
      double acc;
      uint32_t irow, icol, ii;
      arma::rowvec work(d);
            
      uint32_t coreSeed = static_cast<uint32_t>(seeds[0]);
      
      // (Optionally) over-writing the seed here
      #ifdef SUPPORT_OPENMP
      coreSeed = static_cast<uint32_t>( seeds[omp_get_thread_num()] );
      #endif
      
      sitmo::prng_engine engine( coreSeed );
      std::normal_distribution<> normal(0.0, 1.0);
      
      // Filling "out" with standard normal rvs
      #ifdef SUPPORT_OPENMP
      #pragma omp for schedule(static)
      #endif
      for (irow = 0; irow < n; irow++) 
        for (icol = 0; icol < d; icol++) 
           out(irow, icol) = normal(engine);
      
      // Multiplying "out"" by cholesky decomposition of covariance and adding the
      // mean to obtain the desired multivariate normal data rvs.
      #ifdef SUPPORT_OPENMP
      #pragma omp for schedule(static)
      #endif
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
      
      #ifdef SUPPORT_OPENMP
      }
      #endif
      
      return out;
            
    } catch( std::exception& __ex__){
      forward_exception_to_r(__ex__);
    } catch(...){
      ::Rf_error( "c++ exception (unknown reason)" );
    }
    return wrap(NA_REAL);
}