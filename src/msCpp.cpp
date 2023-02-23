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
 *  Mean-shift algorithm
 *
 * See ?ms() for a description of the arguments and output.
*/

RcppExport SEXP msCpp(SEXP init_, SEXP X_, SEXP cholDec_, SEXP ncores_, SEXP tol_, SEXP store_)
{
    using namespace Rcpp;
    
    try{
      
      arma::vec init = as<arma::vec>(init_);
      arma::mat X = as<arma::mat>(X_); 
      arma::mat cholDec = as<arma::mat>(cholDec_); 
      unsigned int ncores = as<unsigned int>(ncores_);
      double tol = as<double>(tol_);
      bool store = as<bool>(store_);
      
      unsigned int d = init.n_elem;
      uint32_t n = X.n_rows;
      
      if(ncores == 0) stop("ncores has to be positive.");
      if( d != X.n_cols ) stop( "The ncol(X) has to equal to length(init)" );
      if( d != cholDec.n_cols ) stop( "The ncol(X) has to equal to cholDec.n_cols" );
      if( d != cholDec.n_rows ) stop( "The ncol(X) has to equal to cholDec.n_rows" );
      
      // Here we set the number of OMP threads, but before we save the original
      // number of threads, so we can re-set before returning.
      int ncores_0;
      #ifdef _OPENMP
      #pragma omp parallel num_threads(1)
      {
       #pragma omp single
       ncores_0 = omp_get_num_threads();
      }
      omp_set_num_threads(ncores);
      #endif
      
      arma::vec currPos = init; 
      arma::vec oldPos(d);
      arma::vec weights(n);
      std::list< std::vector<double> > traj;
            
      if(store) traj.push_back( as< std::vector<double> >( wrap(init) ) );
      
      do{
        
        oldPos = currPos; 
        weights = dmvtInt(X, oldPos, cholDec, false, -1.0, ncores);
        currPos = arma::conv_to<arma::vec>::from( weights.t() * X / sum(weights) ); 
    
        if(store) traj.push_back( as< std::vector<double> >( wrap(currPos) ) );
        
      }
      while( any( abs(currPos - oldPos) > tol ) );
      
      #ifdef _OPENMP
       omp_set_num_threads(ncores_0);
      #endif
            
      return List::create( _["final"] = wrap( currPos ), _["traj"] = wrap( traj ) );
      
    } catch( std::exception& __ex__){
      forward_exception_to_r(__ex__);
    } catch(...){
      ::Rf_error( "c++ exception (unknown reason)" );
    }
    return wrap(NA_REAL);
}