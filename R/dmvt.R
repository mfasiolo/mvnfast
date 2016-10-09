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
#' @param sigma scale matrix (d x d). Alternatively it can be the cholesky decomposition
#'              of the scale matrix. In that case isChol should be set to TRUE. Notice that ff the degrees of 
#'              freedom (the argument \code{df}) is larger than 2, the \code{Cov(X)=sigma*df/(df-2)}.
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
  
  if( !is.matrix(X) ) X <- matrix(X, 1, length(X))
  
  if( !is.matrix(sigma) ) sigma <- as.matrix( sigma )
  
  if(df <= 0.0) stop("df must be positive.")
  
  .Call( "dmvtCpp", 
         X_ = X, 
         mu_ = mu, 
         sigma_ = sigma, 
         df_ = df,
         log_ = log, 
         ncores_ = ncores,
         isChol_ = isChol, 
         PACKAGE = "mvnfast" )
}
