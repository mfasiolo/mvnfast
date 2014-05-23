/* 
 *  Mahalanobis distance functions
 */ 

#ifndef _MAHALANOBIS_H
#define _MAHALANOBIS_H

SEXP mahaCpp(SEXP X, SEXP mu, SEXP sigma, SEXP ncores, SEXP isChol);

arma::vec mahaInt(arma::mat & X,  arma::vec & mu, arma::mat & sigma, int ncores, bool isChol);

arma::vec dmvnInt( arma::mat X, arma::vec mu, arma::mat cholDec, bool log, int ncores);
                  
#endif
