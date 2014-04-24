\name{meanShift}
\alias{meanShift}
\title{Mean-shift mode seeking algorithm}
\usage{
meanShift(X, init, H, tol = 1e-06, traj = FALSE)
}
\arguments{
  \item{X}{n by d matrix containing the data.}

  \item{init}{d-dimensional vector containing the initial
  point for the optimization. By default it is equal to
  \code{colMeans(X)}.}

  \item{H}{Positive definite bandwidth matrix representing
  the covariance of each component of the Gaussian kernel
  density.}

  \item{tol}{Tolerance used to assess the convergence of
  the algorithm, which is stopped if the absolute values of
  increments along all the dimensions are smaller then tol
  at any iteration. Default value is 1e-6.}

  \item{traj}{If \code{FALSE} only the latest iteration is
  returned, if \code{TRUE} the function will return a
  matrix where the i-th row is the position of the
  algorithms at the i-th iteration.}
}
\value{
A list where \code{estim} is a d-dimensional vector
containing the last position of the algorithm, while
\code{allTraj} is a matrix with d-colums representing the
trajectory of the algorithm along each dimension. If
\code{traj == FALSE} the whole trajectory is not stored and
\code{allTraj = NULL}.
}
\description{
Given a sample from a d-dimensional distribution, an
initialization point and a bandwidth the algorithm finds
the nearest mode of the corresponding Gaussian kernel
density.
}
\examples{
# 2 dimensional example
set.seed(434)
X <- matrix(rnorm(400), 100, 2) * c(1, 2)
out <- meanShift(X, init = c(2, 2), H = diag(0.2, 2), traj = TRUE)
plot(X, xlab = "X1", ylab = "X2")
lines(out$allTraj[ , 1], out$allTraj[ , 2], col = 2, lwd = 2)
points(0, 0, col = 3, pch = 3, lwd = 3) # true mode
points(out$estim[1], out$estim[2], col = 4, pch = 3, lwd = 3) # final estimate
}
\author{
Matteo Fasiolo <matteo.fasiolo@gmail.com>.
}

