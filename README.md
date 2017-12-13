
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/mvnfast)](https://cran.r-project.org/package=mvnfast)
[![Build Status](https://travis-ci.org/mfasiolo/mvnfast.svg?branch=master)](https://travis-ci.org/mfasiolo/mvnfast)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/mfasiolo/mvnfast?branch=master&svg=true)](https://ci.appveyor.com/project/mfasiolo/mvnfast)

# **mvnfast**: fast multivariate normal and Student's t distributions

The `mvnfast` R package provides computationally efficient tools related to the multivariate normal and Student's t distributions. The tools are generally faster than those provided by other packages, thanks to the use of C++ code through the `Rcpp`\\`RcppArmadillo` packages and parallelization through the `OpenMP` API. The most important functions are:

- `rmvn()`: simulates multivariate normal random vectors.
- `rmvt()`: simulates Student's t normal random vectors.
- `dmvn()`: evaluates the probability density function of a multivariate normal distribution.  
- `dmvt()`: evaluates the probability density function of a multivariate Student's t distribution.  
- `maha()`: evaluates mahalanobis distances.

See the [vignette](https://mfasiolo.github.io/mvnfast/articles/mvnfast.html) for an introduction to mvnfast and some performance benchmarking.