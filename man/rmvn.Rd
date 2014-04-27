\name{rmvn}
\alias{rmvn}
\title{Fast computation of squared mahalanobis distance}
\usage{
rmvn(n, mu, sigma, ncores = 1, isChol = FALSE)
}
\arguments{
  \item{X}{matrix n by d where each row is a d dimensional
  random vector. Alternatively \code{X} can be a
  d-dimensional vector.}

  \item{mu}{vector of length d, representing the central
  position.}

  \item{sigma}{covariance matrix (d x d). Alternatively is
  can be the cholesky decomposition of the covariance. In
  that case \code{isChol} should be set to \code{TRUE}.}

  \item{isChol}{boolean set to \code{TRUE} is sigma is the
  cholesky decomposition of the covariance.}
}
\value{
a vector of length n where the i-the entry contains the
square mahalanobis distance i-th random vector.
}
\description{
Fast computation of squared mahalanobis distance
}
\examples{
N <- 100
d <- 5
mu <- 1:d
X <- t(t(matrix(rnorm(N*d), N, d)) + mu)
tmp <- matrix(rnorm(d^2), d, d)
mcov <- tcrossprod(tmp, tmp)
myChol <- chol(mcov)

rmvn(10, 1:d, mcov)

\dontrun{
# Performance comparison
library(microbenchmark)

a <- cbind(
  maha(X, mu, mcov),
  maha(X, mu, myChol, isChol = TRUE),
  mahalanobis(X, mu, mcov))

# Same output as mahalanobis
a[ , 1] / a[, 3]
a[ , 2] / a[, 3]

microbenchmark(maha(X, mu, mcov),
               maha(X, mu, myChol, isChol = TRUE),
               mahalanobis(X, mu, mcov))
}
}
\author{
Matteo Fasiolo <matteo.fasiolo@gmail.com>
}

