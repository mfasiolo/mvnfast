/* 
 *  Mahalanobis distance functions
 */ 

#ifndef _MAHALANOBIS_H
#define _MAHALANOBIS_H

SEXP mahaCpp(SEXP X, SEXP mu, SEXP sigma, SEXP isChol);
arma::rowvec mahaInt(arma::mat & X,  arma::vec & mu, arma::mat & sigma, bool isChol);
                  
#endif
