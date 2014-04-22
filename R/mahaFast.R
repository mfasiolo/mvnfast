######
## Fast computation of mahalanobis distance
######
#' Fast computation of squared mahalanobis distance
#'
#' @param X matrix n by d where each row is a d dimensional random vector. Alternatively \code{X} can be a d-dimensional vector.
#' @param mu vector of length d, representing the central position.
#' @param sigma covariance matrix (d x d). Alternatively is can be the cholesky decomposition
#'              of the covariance. In that case \code{isChol} should be set to \code{TRUE}.
#' @param isChol boolean set to \code{TRUE} is sigma is the cholesky decomposition of the covariance.
#' @return a vector of length n where the i-the entry contains the square mahalanobis distance i-th random vector.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>
#' @examples
#' N <- 100
#' d <- 5
#' mu <- 1:d
#' X <- t(t(matrix(rnorm(N*d), N, d)) + mu)
#' tmp <- matrix(rnorm(d^2), d, d)
#' mcov <- tcrossprod(tmp, tmp)
#' myChol <- chol(mcov)
#' 
#' rbind(head(mahaFast(X, mu, mcov), 10),
#'       head(mahaFast(X, mu, myChol, isChol = TRUE), 10),
#'       head(mahalanobis(X, mu, mcov), 10))
#' 
#' \dontrun{
#' # Performance comparison
#' library(microbenchmark)
#' 
#' a <- cbind(
#'   mahaFast(X, mu, mcov),
#'   mahaFast(X, mu, myChol, isChol = TRUE),
#'   mahalanobis(X, mu, mcov))
#'   
#' # Same output as mahalanobis
#' a[ , 1] / a[, 3]
#' a[ , 2] / a[, 3]
#' 
#' microbenchmark(mahaFast(X, mu, mcov),
#'                mahaFast(X, mu, myChol, isChol = TRUE),
#'                mahalanobis(X, mu, mcov))
#' }
#' @export 

mahaFast <- function(X, mu, sigma, isChol = FALSE)
{
  if( !is.matrix(X) ) X <- matrix(X, 1, length(X))
    
  drop(.Call( "mahaCpp", 
              X_ = X, 
              mu_ = mu, 
              sigma_ = sigma, 
              isChol_ = isChol, 
              PACKAGE = "synlik" )) 
}


