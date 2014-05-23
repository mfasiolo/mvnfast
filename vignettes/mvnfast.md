
<!--
%\VignetteEngine{knitr::docco_linear}
%\VignetteIndexEntry{synlik_vignette}
-->
  
An Introduction to **mvnfast**
=======================================
  




Introduction
------------
  
The `mvn` R package provides computationally efficient tools related to the multivariate normal distribution. 
The tools are generally faster than those provided by other packages, thanks to the use of C++ code through the 
`Rcpp`\\`RcppArmadillo` packages and parallelization through the `OpenMP` API. The most important functions are:
- `rmvn()`: simulates multivariate normal random vectors.
- `dmvn()`: evaluates the probability density function of a multivariate normal distribution.  
- `maha()`: evaluates mahalanobis distances.

In the following sections we will benchmark each function against equivalent functions provided by other packages, while in the final section we provide an example application.


Simulating multivariate normal random vectors
----------------------------
  
Simulating multivariate normal random variables is an essential step in many Monte Carlo algorithms (such as MCMC or Particle Filters),
 hence this operations has to be as fast as possible. Here we compare the `rmvn` function with the equivalent function `rmvnorm` 
(from the `mvtnorm` package) and `mvrnorm` (from the `MASS` package). In particular, we simulate $10^4$ twenty-dimensional random vectors:

```r
library("microbenchmark")
library("mvtnorm")
library("mvnfast")
library("MASS")

N <- 10000
d <- 20

# Creating mean and covariance matrix
mu <- 1:d
tmp <- matrix(rnorm(d^2), d, d)
mcov <- tcrossprod(tmp, tmp)

microbenchmark(rmvn(N, mu, mcov, ncores = 2),
               rmvn(N, mu, mcov),
               rmvnorm(N, mu, mcov),
               mvrnorm(N, mu, mcov))
```

```
## Unit: milliseconds
##                           expr   min    lq median    uq    max neval
##  rmvn(N, mu, mcov, ncores = 2) 19.20 19.44  20.26 21.43  45.25   100
##              rmvn(N, mu, mcov) 37.39 37.57  38.35 39.30  63.65   100
##           rmvnorm(N, mu, mcov) 48.71 55.62  59.51 66.93  93.03   100
##           mvrnorm(N, mu, mcov) 46.47 52.24  55.66 61.40 108.49   100
```


In this example `rmvn` cuts the computational time, relative to the alternatives, even when a single core is used. This gain is attributable to several factors: the use of C++ code and efficient numerical algorithms to simulate the random variables. Parallelizing the computation on two cores gives another appreciable speed-up. To be fair, it is necessary to point out that `rmvnorm` and `mvrnorm` have many more safety check on the user's input than `rmvn`. This is true also for the functions described in the next sections.

Finally, notice that this function does not use one of the Random Number Generators (RNGs) provided by R, but one 
of the parallel cryptographic RNGs described in (Salmon et al., 2011) and available [here](http://www.sitmo.com/article/parallel-random-number-generator-in-c/). It is important to point out that this RNG can safely be used in parallel, without risk of collisions between parallel sequence of random numbers, as detailed in the above reference.

Evaluating the multivariate normal density
----------------------------

Here we compare the `dmvn` function, which evaluates the multivariate normal density,  with the equivalent function `dmvtnorm` (from the `mvtnorm` package). 
In particular we evaluate the density of $10^4$ twenty-dimensional random vectors:

```r
# Generating random vectors 
N <- 10000
d <- 20
mu <- 1:d
tmp <- matrix(rnorm(d^2), d, d)
mcov <- tcrossprod(tmp, tmp)
X <- rmvn(N, mu, mcov)

microbenchmark(dmvn(X, mu, mcov, ncores = 2),
               dmvn(X, mu, mcov),
               dmvnorm(X, mu, mcov, trustme = T))
```

```
## Unit: milliseconds
##                               expr   min    lq median    uq    max neval
##      dmvn(X, mu, mcov, ncores = 2) 5.114 5.203  5.234 5.282  6.570   100
##                  dmvn(X, mu, mcov) 7.256 7.318  7.349 7.404  7.838   100
##  dmvnorm(X, mu, mcov, trustme = T) 7.688 8.835  9.020 9.186 33.993   100
```

Again, we get some speed-up using C++ code and some more from the parallelization.

Evaluating the Mahalanobis distance
----------------------------

Finally, we compare the `maha` function, which evaluates the square [mahalanobis distance](http://en.wikipedia.org/wiki/Mahalanobis_distance) with the equivalent function `mahalanobis` (from the `stats` package). 
Also in the case we use $10^4$ twenty-dimensional random vectors:

```r
# Generating random vectors 
N <- 10000
d <- 20
mu <- 1:d
tmp <- matrix(rnorm(d^2), d, d)
mcov <- tcrossprod(tmp, tmp)
X <- rmvn(N, mu, mcov)

microbenchmark(maha(X, mu, mcov, ncores = 2),
               maha(X, mu, mcov),
               mahalanobis(X, mu, mcov))
```

```
## Unit: milliseconds
##                           expr    min     lq median     uq    max neval
##  maha(X, mu, mcov, ncores = 2)  3.368  3.467  3.774  4.741  8.225   100
##              maha(X, mu, mcov)  5.535  5.593  5.771  7.048 31.478   100
##       mahalanobis(X, mu, mcov) 11.892 15.436 18.157 21.886 53.244   100
```

The acceleration is similar to that obtained in the previous sections.

Example: mean-shift mode seeking algorithm
----------------------------

As an example application of the `dmvn` function, we implemented the [mean-shift mode seeking](http://en.wikipedia.org/wiki/Mean-shift) algorithm.
This procedure can be used to find the mode or maxima of a kernel density function, and it can be used to set up
clustering algorithms. Here we simulate $10^4$ d-dimensional random vectors from mixture of normal distributions: 

```r
set.seed(5135)
N <- 10000
d <- 2
mu1 <- c(0, 0); mu2 <- c(2, 3)
Cov1 <- matrix(c(1, 0, 0, 2), 2, 2)
Cov2 <- matrix(c(1, -0.9, -0.9, 1), 2, 2)

bin <- rbinom(N, 1, 0.5)

X <- bin * rmvn(N, mu1, Cov1) + (!bin) * rmvn(N, mu2, Cov2)
```

Finally, we plot the resulting probability density and, starting from 10 initial points,  we use mean-shift to converge to the nearest mode:

```r
# Plotting
np <- 100
xvals <- seq(min(X[ , 1]), max(X[ , 1]), length.out = np)
yvals <- seq(min(X[ , 2]), max(X[ , 2]), length.out = np)
theGrid <- expand.grid(xvals, yvals) 
theGrid <- as.matrix(theGrid)
dens <- 0.5 * dmvn(theGrid, mu1, Cov1) + 0.5 * dmvn(theGrid, mu2, Cov2)
plot(X[ , 1], X[ , 2], pch = '.', lwd = 0.01, col = 3)
contour(x = xvals, y = yvals, z = matrix(dens, np, np),
        levels = c(0.002, 0.01, 0.02, 0.04, 0.08, 0.15 ), add = TRUE, lwd = 2)

# Mean-shift
library(plyr)
inits <- matrix(c(-2, 2, 0, 3, 4, 3, 2, 5, 2, -3, 2, 2, 0, 2, 3, 0, 0, -4, -2, 6), 
                10, 2, byrow = TRUE)
traj <- alply(inits,
              1,
              function(input)
                  ms(X = X, 
                     init = input, 
                     H = 0.05 * cov(X), 
                     ncores = 2, 
                     store = TRUE)$traj
              )

invisible( lapply(traj, 
                  function(input){ 
                    lines(input[ , 1], input[ , 2], col = 2, lwd = 1.5)
                    points(tail(input[ , 1]), tail(input[ , 2]))
           }))
```

<img src="figure/mixPlot.png" title="plot of chunk mixPlot" alt="plot of chunk mixPlot" style="display:block; margin: auto" style="display: block; margin: auto;" />

As we can see from the plot, each initial point leads one of two points that are very close to the true mode. Notice that the bandwidth for the kernel density estimator was chosen by trial-and-error, and less arbitrary choices are certainly possible in real applications. 
 
References
----------------------------
  
  * Dirk Eddelbuettel and Romain Francois (2011). Rcpp: Seamless R and C++ Integration. Journal of Statistical Software, 40(8),
  1-18. URL http://www.jstatsoft.org/v40/i08/.
  
  * Eddelbuettel, Dirk (2013) Seamless R and C++ Integration with Rcpp. Springer, New York. ISBN 978-1-4614-6867-7.
  
  *  Dirk Eddelbuettel, Conrad Sanderson (2014). RcppArmadillo: Accelerating R with high-performance C++ linear algebra. Computational
  Statistics and Data Analysis, Volume 71, March 2014, pages 1054-1063. URL http://dx.doi.org/10.1016/j.csda.2013.02.005

  * http://openmp.org/
  
  * John K. Salmon, Mark A. Moraes, Ron O. Dror, and David E. Shaw (2011). Parallel Random Numbers: As Easy as 1, 2, 3.
    D. E. Shaw Research, New York, NY 10036, USA.
















