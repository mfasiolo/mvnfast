#########
#' Fast computation of the multivariate Student's t density.
#'
#' @param X matrix n by d where each row is a d dimensional random vector. Alternatively \code{X} can be a d-dimensional vector.
#' @param mu vector of length d, representing the mean of the distribution.
#' @param sigma scale matrix (d x d). Alternatively it can be the cholesky decomposition
#'              of the scale matrix. In that case isChol should be set to TRUE. Notice that ff the degrees of 
#'              freedom (the argument \code{df}) is larger than 2, the \code{Cov(X)=sigma*df/(df-2)}.
#' @param df a positive scalar representing the degrees of freedom.
#' @param log boolean set to true the logarithm of the pdf is required.
#' @param ncores Number of cores used. The parallelization will take place only if OpenMP is supported.
#' @param isChol boolean set to true is \code{sigma} is the cholesky decomposition of the covariance matrix.
#' @return A vector of length n where the i-the entry contains the pdf of the i-th random vector.
#' @details There are many candidates for the multivariate generalization of Student's t-distribution, here we use
#'          the parametrization described here \url{https://en.wikipedia.org/wiki/Multivariate_t-distribution}. NB: at the moment 
#'          the parallelization does not work properly on Solaris OS when \code{ncores>1}. Hence, \code{dmvt()} checks if the OS 
#'          is Solaris and, if this the case, it imposes \code{ncores==1}. 
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> 
#' @examples
#' N <- 100
#' d <- 5
#' mu <- 1:d
#' df <- 4
#' X <- t(t(matrix(rnorm(N*d), N, d)) + mu)
#' tmp <- matrix(rnorm(d^2), d, d)
#' mcov <- tcrossprod(tmp, tmp)  + diag(0.5, d)
#' myChol <- chol(mcov)
#' 
#' head(dmvt(X, mu, mcov, df = df), 10)
#' head(dmvt(X, mu, myChol, df = df, isChol = TRUE), 10)
#' 
#' @export dmvt

dmvt <- function(X, mu, sigma, df, log = FALSE, ncores = 1, isChol = FALSE){
  
  if( !is.numeric(mu) ) mu <- as.numeric(mu)
  
  if( !is.matrix(X) ) X <- matrix(X, 1, length(X))
  
  if( !is.matrix(sigma) ) sigma <- as.matrix( sigma )
  
  if(df <= 0.0) stop("df must be positive.")
  
  if(df == Inf){ df <- -1 } # Fall back on multivariate normal case
  
  if( ncores > 1 && grepl('SunOS', Sys.info()['sysname']) ){
    
    message("dmvt() cannot be used on multiple cores under Solaris. I am resetting \"ncores\" to 1.")
    
    ncores <- 1
    
  }
  
  .Call( "dmvtCpp", 
         X_ = X, 
         mu_ = mu, 
         sigma_ = sigma, 
         df_ = df,
         log_ = log, 
         ncores_ = ncores,
         isChol_ = isChol)
}
