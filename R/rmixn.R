##############################################################
#' Fast simulation of r.v.s from a mixture of multivariate normal densities
#'
#' @param n number of random vectors to be simulated.
#' @param mu an (m x d) matrix, where m is the number of mixture components.
#' @param sigma as list of m covariance matrices (d x d) on for each mixture component. 
#'              Alternatively it can be a list of m cholesky decomposition of the covariance. 
#'              In that case \code{isChol} should be set to \code{TRUE}.
#' @param w vector of length m, containing the weights of the mixture components.
#' @param ncores Number of cores used. The parallelization will take place only if OpenMP is supported.
#' @param isChol boolean set to true is \code{sigma} is the cholesky decomposition of the covariance matrix.
#' @param retInd when set to \code{TRUE} an attribute called "index" will be added to the output matrix of random variables.
#'               This is a vector specifying to which mixture components each random vector belongs. \code{FALSE} by default.
#' @param A an (optional) numeric matrix of dimension (n x d), which will be used to store the output random variables.
#'        It is useful when n and d are large and one wants to call \code{rmvn()} several times, without reallocating memory
#'        for the whole matrix each time. NB: the element of \code{A} must be of class "numeric".
#' @param kpnames if \code{TRUE} the dimensions' names are preserved. That is, the i-th column of the output
#'                has the same name as the i-th entry of \code{mu} or the i-th column of \code{sigma}. 
#'                \code{kpnames==FALSE} by default.
#' @return If \code{A==NULL} (default) the output is an (n x d) matrix where the i-th row is the i-th simulated vector.
#'         If \code{A!=NULL} then the random vector are store in \code{A}, which is provided by the user, and the function
#'         returns \code{NULL}. Notice that if \code{retInd==TRUE} an attribute called "index" will be added to A.
#'         This is a vector specifying to which mixture components each random vector belongs.
#' @details Notice that this function does not use one of the Random Number Generators (RNGs) provided by R, but one 
#'          of the parallel cryptographic RNGs described in (Salmon et al., 2011). It is important to point out that this
#'          RNG can safely be used in parallel, without risk of collisions between parallel sequence of random numbers.
#'          The initialization of the RNG depends on R's seed, hence the \code{set.seed()} function can be used to 
#'          obtain reproducible results. Notice though that changing \code{ncores} causes most of the generated numbers
#'          to be different even if R's seed is the same (see example below). NB: at the moment the RNG does not work
#'          properly on Solaris OS when \code{ncores>1}. Hence, \code{rmixn()} checks if the OS is Solaris and, if this the case, 
#'          it imposes \code{ncores==1}. 
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>, C++ RNG engine by Thijs van den Berg <http://sitmo.com/>.
#' @references  John K. Salmon, Mark A. Moraes, Ron O. Dror, and David E. Shaw (2011). Parallel Random Numbers: As Easy as 1, 2, 3.
#'              D. E. Shaw Research, New York, NY 10036, USA.
#' @examples
#' # Create mixture of two components
#' mu <- matrix(rep(c(1, 2, 10, 20), 2), 2, 2, byrow = TRUE)
#' sigma <- list(diag(c(1, 10)), matrix(c(1, -0.9, -0.9, 1), 2, 2))
#' w <- c(0.1, 0.9)
#' 
#' # Simulate
#' X <- rmixn(1e4, mu, sigma, w, retInd = TRUE)
#' plot(X, pch = '.', col = attr(X, "index"))
#' 
#' # Simulate with fixed seed
#' set.seed(414)
#' rmixn(4, mu, sigma, w)
#' 
#' set.seed(414)
#' rmixn(4, mu, sigma, w)
#' 
#' set.seed(414)  
#' rmixn(4, mu, sigma, w, ncores = 2) # r.v. generated on the second core are different
#' 
#' ###### Here we create the matrix that will hold the simulated random variables upfront.
#' A <- matrix(NA, 4, 2)
#' class(A) <- "numeric" # This is important. We need the elements of A to be of class "numeric". 
#' 
#' set.seed(414)
#' rmixn(4, mu, sigma, w, ncores = 2, A = A) # This returns NULL ...
#' A                                         # ... but the result is here
#' 
#' @export rmixn
#'
rmixn <- function(n, mu, sigma, w, ncores = 1, isChol = FALSE, retInd = FALSE, A = NULL, kpnames = FALSE)
{
  d <- ncol(mu)
  m <- length(w)
  
  if( !is.matrix(sigma) ) sigma <- as.matrix( sigma )
  
  if( ncores > 1 && grepl('SunOS', Sys.info()['sysname']) ){
    
    message("rmvn() cannot be used on multiple cores under Solaris. I am resetting \"ncores\" to 1.")
    
    ncores <- 1
    
  }
  
  # Create output matrix
  if( is.null(A) ) {
    retMat <- TRUE # We return a matrix
    A <- matrix(nrow = n, ncol = d) 
    class(A) <- "numeric"
  } else {
    retMat <- FALSE # We return NULL
    if( class(A[1, 1]) != "numeric" ){ 
      stop("class(A[1, 1]) != \"numeric\", to avoid this do class(A)<-\"numeric\".")
    }
  } 
  
  # Associate each sample with a mixture component
  indV <- sample(0:(m-1), n, prob = w, replace = T)

  .Call( "rmixnCpp", 
         n_ = n, 
         mu_ = mu, 
         sigma_ = sigma,
         indV_ = indV,
         ncores_ = ncores,
         isChol_ = isChol, 
         retInd_ = retInd,
         A_ = A )
  
  if( kpnames ) # (Optionally add dimensions names)
  {
    if( is.null(nam <- dimnames(mu)[[2L]]) ){ nam <- dimnames(sigma[[1L]])[[1L]] }
    if( !is.null(nam) ){ dimnames(A) <- list(NULL, nam) }
  }
  
  # Return a matrix if no storage was provided and NULL if it was provided.
  if( retMat ) {
    return( A );
  } else {
    return( invisible(NULL) )
  }
  
}