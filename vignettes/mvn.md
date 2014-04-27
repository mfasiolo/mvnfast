
<!--
%\VignetteEngine{knitr::docco_linear}
%\VignetteIndexEntry{synlik_vignette}
-->
  
An Introduction to **mvn**
=======================================
  




Introduction
------------
  
The `mvn` R package provides computationally efficient tools related to the multivariate normal distribution. 
The tools are generally faster than those provided by other packages, thanks to the use of C++ code through the 
`Rcpp` package and of parallel computation through the `OpenMP` API.


Evaluating the multivariate normal density
----------------------------

Computing the density of multivariate normal is an essential step of MCMC or clustering algorithms, hence this
operations has to be as fast as possible. Here we compare the `dmvn` function with the equivalent function `dmvtnorm`, 
provided by the `mvtnorm` package. We start by evaluating the density of $10^4$


```r
library("microbenchmark")
library("mvtnorm")
library("mvn")

N <- 10000
d <- 50
mu <- 1:d
X <- t(t(matrix(rnorm(N*d), N, d)) + mu)
tmp <- matrix(rnorm(d^2), d, d)
mcov <- tcrossprod(tmp, tmp)
myChol <- chol(mcov)

```




References
----------------------------
  
  * Simon N Wood. Statistical inference for noisy nonlinear ecological dynamic systems. Nature, 466(7310):1102--1104, 2010.



















