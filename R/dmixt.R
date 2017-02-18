##############################################################
#' Fast density computation for mixture of multivariate Student's t distributions.
#'
#' @param X matrix n by d where each row is a d dimensional random vector. Alternatively \code{X} can be a d-dimensional vector.
#' @param mu an (m x d) matrix, where m is the number of mixture components.
#' @param sigma as list of m covariance matrices (d x d) on for each mixture component. 
#'              Alternatively it can be a list of m cholesky decomposition of the covariance. 
#'              In that case \code{isChol} should be set to \code{TRUE}.
#' @param df a positive scalar representing the degrees of freedom. All the densities in the mixture have the same \code{df}.
#' @param w vector of length m, containing the weights of the mixture components.
#' @param log boolean set to true the logarithm of the pdf is required.
#' @param ncores Number of cores used. The parallelization will take place only if OpenMP is supported.
#' @param isChol boolean set to true is \code{sigma} is the cholesky decomposition of the covariance matrix.
#' @param A an (optional) numeric matrix of dimension (m x d), which will be used to store the evaluations of each mixture
#'        density over each mixture component. It is useful when m and n are large and one wants to call \code{dmixt()} 
#'        several times, without reallocating memory for the whole matrix each time. NB1: \code{A} will be modified, 
#'        not copied! NB2: the element of \code{A} must be of class "numeric".
#' @return A vector of length n where the i-the entry contains the pdf of the i-th random vector (i.e. the i-th row of \code{X}).
#' @details There are many candidates for the multivariate generalization of Student's t-distribution, here we use
#'          the parametrization described here \url{https://en.wikipedia.org/wiki/Multivariate_t-distribution}. NB: at the moment 
#'          the parallelization does not work properly on Solaris OS when \code{ncores>1}. Hence, \code{dmixt()} checks if the OS 
#'          is Solaris and, if this the case, it imposes \code{ncores==1}. 
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#' @examples
#' #### 1) Example use
#' # Set up mixture density
#' df <- 6
#' mu <- matrix(rep(c(1, 2, 10, 20), 2), 2, 2, byrow = TRUE)
#' sigma <- list(diag(c(1, 10)), matrix(c(1, -0.9, -0.9, 1), 2, 2))
#' w <- c(0.1, 0.9)
#' 
#' # Simulate
#' X <- rmixt(1e4, mu, sigma, df, w)
#' 
#' # Evaluate density
#' ds <- dmixt(X, mu, sigma, w = w, df = df)
#' head(ds)
#' 
#' ##### 2) More complicated example
#' # Define mixture
#' set.seed(5135)
#' N <- 10000
#' d <- 2
#' df = 10
#' w <- rep(1, 2) / 2
#' mu <- matrix(rep(c(0, 0, 2, 3), 2), 2, 2, byrow = TRUE) 
#' sigma <- list(matrix(c(1, 0, 0, 2), 2, 2), matrix(c(1, -0.9, -0.9, 1), 2, 2)) 
#' 
#' # Simulate random variables
#' X <- rmixt(N, mu, sigma, w = w, df = df, retInd = TRUE)
#' 
#' # Plot mixture density
#' np <- 100
#' xvals <- seq(min(X[ , 1]), max(X[ , 1]), length.out = np)
#' yvals <- seq(min(X[ , 2]), max(X[ , 2]), length.out = np)
#' theGrid <- expand.grid(xvals, yvals) 
#' theGrid <- as.matrix(theGrid)
#' dens <- dmixt(theGrid, mu, sigma, w = w, df = df)
#' plot(X, pch = '.', col = attr(X, "index")+1)
#' contour(x = xvals, y = yvals, z = matrix(dens, np, np),
#'         levels = c(0.002, 0.01, 0.02, 0.04, 0.08, 0.15 ), add = TRUE, lwd = 2)
#' 
#' @export dmixt
#'
dmixt <- function(X, mu, sigma, df, w, log = FALSE, ncores = 1, isChol = FALSE, A = NULL)
{
  d <- ncol(mu)
  m <- length(w)
  n <- nrow(X)
  
  if( !is.matrix(sigma) ) sigma <- as.matrix( sigma )
  
  if(df <= 0.0) stop("df must be positive.")
  
  if(df == Inf){ df <- -1 } # Fall back on multivariate normal case 
  
  if( ncores > 1 && grepl('SunOS', Sys.info()['sysname']) ){
    
    message("rmvn() cannot be used on multiple cores under Solaris. I am resetting \"ncores\" to 1.")
    
    ncores <- 1
    
  }
  
  # Create output matrix
  if( is.null(A) ) {
    A <- matrix(nrow = n, ncol = m) 
    class(A) <- "numeric"
  } else {
    if( class(A[1, 1]) != "numeric" ){ 
      stop("class(A[1, 1]) != \"numeric\", to avoid this do class(A)<-\"numeric\".")
    }
  } 
  
  out <-.Call( "dmixtCpp", 
               X_ = X, 
               mu_ = mu, 
               sigma_ = sigma,
               df_ = df,
               w_ = w,
               log_ = log,
               ncores_ = ncores,
               isChol_ = isChol, 
               A_ = A )
  
  return( out )
  
}