
test_that("Checking dmvn against dmvnorm", {
  library("mvtnorm")
  
  ###### d = 1, n = 1 case
  N <- 1
  d <- 1
  mu <- 1:d
  X <- t(t(matrix(rnorm(N*d), N, d)) + mu)
  tmp <- matrix(rnorm(d^2), d, d)
  mcov <- tcrossprod(tmp, tmp)
  myChol <- chol(mcov)
  
  a <- cbind(
    dmvn(X, mu, mcov),
    dmvn(X, mu, myChol, isChol = TRUE),
    dmvnorm(X, mu, mcov))
  
  expect_less_than(abs(a[ , 1] - a[, 3]), 1e-10)
  expect_less_than(abs(a[ , 2] - a[, 3]), 1e-10)
  
  ###### d = 1, n = 100 case
  N <- 100
  d <- 1
  mu <- 1:d
  X <- t(t(matrix(rnorm(N*d), N, d)) + mu)
  tmp <- matrix(rnorm(d^2), d, d)
  mcov <- tcrossprod(tmp, tmp)
  myChol <- chol(mcov)
  
  a <- cbind(
    dmvn(X, mu, mcov),
    dmvn(X, mu, myChol, isChol = TRUE),
    dmvnorm(X, mu, mcov))
  
  expect_less_than(sum(abs(a[ , 1] - a[, 3])), 1e-10)
  expect_less_than(sum(abs(a[ , 2] - a[, 3])), 1e-10)
  
  ###### d = 10, n = 1 case
  N <- 1
  d <- 10
  mu <- 1:d
  X <- t(t(matrix(rnorm(N*d), N, d)) + mu)
  tmp <- matrix(rnorm(d^2), d, d)
  mcov <- tcrossprod(tmp, tmp)
  myChol <- chol(mcov)
  
  a <- cbind(
    dmvn(X, mu, mcov),
    dmvn(X, mu, myChol, isChol = TRUE),
    dmvnorm(X, mu, mcov))
  
  expect_less_than(sum(abs(a[ , 1] - a[, 3])), 1e-10)
  expect_less_than(sum(abs(a[ , 2] - a[, 3])), 1e-10)
  
  ###### d = 10, n = 100 case
  N <- 100
  d <- 10
  mu <- 1:d
  X <- t(t(matrix(rnorm(N*d), N, d)) + mu)
  tmp <- matrix(rnorm(d^2), d, d)
  mcov <- tcrossprod(tmp, tmp)
  myChol <- chol(mcov)
  
  a <- cbind(
    dmvn(X, mu, mcov),
    dmvn(X, mu, myChol, isChol = TRUE),
    dmvnorm(X, mu, mcov))
  
  expect_less_than(sum(abs(a[ , 1] - a[, 3])), 1e-10)
  expect_less_than(sum(abs(a[ , 2] - a[, 3])), 1e-10)
  
})