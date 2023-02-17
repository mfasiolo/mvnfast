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
 * Simulate random variables from a multivariate normal distribution
 *
 * See ?rmvn() for a description of the arguments and output.
*/

RcppExport SEXP rmvnCpp(SEXP n_,  
                        SEXP mu_,  
                        SEXP sigma_,
                        SEXP ncores_,
                        SEXP isChol_, 
                        SEXP A_) 
{ 
    using namespace Rcpp;
    
    try{
      
      uint32_t n = as<uint32_t>(n_);
      arma::rowvec mu = as<arma::rowvec>(mu_);  
      arma::mat sigma = as<arma::mat>(sigma_); 
      unsigned int  ncores = as<unsigned int>(ncores_); 
      bool isChol = as<bool>(isChol_); 
      NumericMatrix A = NumericMatrix(A_);
      
      unsigned int d = mu.n_elem;
      
      if(n < 1) stop("n should be a positive integer");
      if(ncores == 0) stop("ncores has to be positive");
      if(d != sigma.n_cols) stop("mu.n_elem != sigma.n_cols");
      if(d != sigma.n_rows) stop("mu.n_elem != sigma.n_rows");
      if(d != A.ncol()) stop("mu.n_elem != A.ncol()");
      if(n != A.nrow()) stop("n != A.nrow()");
      
      #ifdef _OPENMP
       omp_set_num_threads(ncores);
      #endif
      
      // The A matrix that will be filled with firstly with standard normal rvs,
      // and finally with multivariate normal rvs.
      // We A wrap into a arma::mat "tmp" without making a copy.
      arma::mat tmp( A.begin(), A.nrow(), A.ncol(), false );
      
      RNGScope scope; // Declare RNGScope after the output in order to avoid a known Rcpp bug.
    
      // Calculate cholesky dec unless sigma is already a cholesky dec.
      arma::mat cholDec;
      if( isChol == false ) {
        cholDec = trimatu( arma::chol(sigma) );
      }
      else{
        cholDec = trimatu( sigma );
      }
                  
      // What I do to seed the sitmo::prng_engine is tricky. I produce "ncores" uniform numbers between 1 and the largest uint32_t,
      // which are the seeds. I put the first one in "coreSeed". If there is no support for OpenMP only this seed
      // will be used, as the computations will be sequential. If there is support for OpenMP, "coreSeed" will
      // be over-written, so that each core will get its own seed.
      //NumericVector seeds = runif(ncores, 1.0, std::numeric_limits<uint32_t>::max());
      
      NumericVector seeds(ncores);
      seeds[0] = runif(1, 1.0, std::numeric_limits<uint32_t>::max())[0];
      for (unsigned int j = 0;  j < ncores - 1; j++){
        seeds[j+1] = seeds[j] - 1.0;
        if (seeds[j+1] < 1.0) seeds[j] = std::numeric_limits<uint32_t>::max() - 1.0;
      }
      
      #ifdef _OPENMP
      #pragma omp parallel num_threads(ncores) if(ncores > 1)
      {
      #endif
      
      double acc;
      uint32_t irow, icol, ii;
      arma::rowvec work(d);
            
      uint32_t coreSeed = static_cast<uint32_t>(seeds[0]);
      
      // (Optionally) over-writing the seed here
      #ifdef _OPENMP
      coreSeed = static_cast<uint32_t>( seeds[omp_get_thread_num()] );
      #endif
       
      sitmo::prng_engine engine( coreSeed );
      boost::normal_distribution<> normal(0.0, 1.0);
      
      // Filling "out" with standard normal rvs
      #ifdef _OPENMP
      #pragma omp for schedule(static)
      #endif
      for (irow = 0; irow < n; irow++) 
        for (icol = 0; icol < d; icol++) 
           A(irow, icol) = normal(engine);
      
      // Multiplying "out"" by cholesky decomposition of covariance and adding the
      // mean to obtain the desired multivariate normal data rvs.
      #ifdef _OPENMP
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
      
      #ifdef _OPENMP
      }
      #endif
      
      return R_NilValue;
            
    } catch( std::exception& __ex__){
      forward_exception_to_r(__ex__);
    } catch(...){
      ::Rf_error( "c++ exception (unknown reason)" );
    }
    return wrap(NA_REAL);
}