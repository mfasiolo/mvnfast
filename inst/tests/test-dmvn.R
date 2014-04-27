
test_that("Checking dmvn against dmvnorm", {
  library("mvtnorm")
  
  ##########
  ###### d = 1, n = 1 case
  ##########
  N <- 1
  d <- 1
  mu <- 1:d
  X <- t(t(matrix(rnorm(N*d), N, d)) + mu)
  tmp <- matrix(rnorm(d^2), d, d)
  mcov <- tcrossprod(tmp, tmp)
  myChol <- chol(mcov)
  
  ##### dmvn()
  bench <- dmvnorm(X, mu, mcov, log = T)
  # Sequential
  expect_less_than(abs(dmvn(X, mu, mcov, log = T) - bench), 1e-10)
  expect_less_than(abs(dmvn(X, mu, myChol, isChol = TRUE, log = T) - bench), 1e-10)
  # Parallel
  expect_less_than(abs(dmvn(X, mu, mcov, ncores = 2, log = T) - bench), 1e-10)
  expect_less_than(abs(dmvn(X, mu, myChol, ncores = 2, isChol = TRUE, log = T) - bench), 1e-10)
  
  ##### maha()
  bench <- mahalanobis(X, mu, mcov)
  # Sequential
  expect_less_than(abs(maha(X, mu, mcov) - bench), 1e-10)
  expect_less_than(abs(maha(X, mu, myChol, isChol = TRUE) - bench), 1e-10)
  # Parallel
  expect_less_than(abs(maha(X, mu, mcov, ncores = 2) - bench), 1e-10)
  expect_less_than(abs(maha(X, mu, myChol, ncores = 2, isChol = TRUE) - bench), 1e-10)
  
  #############
  ###### d = 1, n = 100 case
  #############
  N <- 100
  d <- 1
  mu <- 1:d
  X <- t(t(matrix(rnorm(N*d), N, d)) + mu)
  tmp <- matrix(rnorm(d^2), d, d)
  mcov <- tcrossprod(tmp, tmp)
  myChol <- chol(mcov)
  
  ##### dmvn()
  bench <- dmvnorm(X, mu, mcov, log = T)
  # Sequential
  expect_less_than(sum(abs(dmvn(X, mu, mcov, log = T) - bench)), 1e-10)
  expect_less_than(sum(abs(dmvn(X, mu, myChol, isChol = TRUE, log = T) - bench)), 1e-10)
  # Parallel
  expect_less_than(sum(abs(dmvn(X, mu, mcov, ncores = 2, log = T) - bench)), 1e-10)
  expect_less_than(sum(abs(dmvn(X, mu, myChol, ncores = 2, isChol = TRUE, log = T) - bench)), 1e-10)
  
  ##### maha()
  bench <- mahalanobis(X, mu, mcov)
  # Sequential
  expect_less_than(sum(abs(maha(X, mu, mcov) - bench)), 1e-10)
  expect_less_than(sum(abs(maha(X, mu, myChol, isChol = TRUE) - bench)), 1e-10)
  # Parallel
  expect_less_than(sum(abs(maha(X, mu, mcov, ncores = 2) - bench)), 1e-10)
  expect_less_than(sum(abs(maha(X, mu, myChol, ncores = 2, isChol = TRUE) - bench)), 1e-10)
  
  #############
  ###### d = 10, n = 10 case
  #############
  N <- 1
  d <- 10
  mu <- 1:d
  X <- t(t(matrix(rnorm(N*d), N, d)) + mu)
  tmp <- matrix(rnorm(d^2), d, d)
  mcov <- tcrossprod(tmp, tmp)
  myChol <- chol(mcov)
  
  ##### dmvn()
  bench <- dmvnorm(X, mu, mcov, log = T)
  # Sequential
  expect_less_than(sum(abs(dmvn(X, mu, mcov, log = T) - bench)), 1e-10)
  expect_less_than(sum(abs(dmvn(X, mu, myChol, isChol = TRUE, log = T) - bench)), 1e-10)
  # Parallel
  expect_less_than(sum(abs(dmvn(X, mu, mcov, ncores = 2, log = T) - bench)), 1e-10)
  expect_less_than(sum(abs(dmvn(X, mu, myChol, ncores = 2, isChol = TRUE, log = T) - bench)), 1e-10)
  
  ##### maha()
  bench <- mahalanobis(X, mu, mcov)
  # Sequential
  expect_less_than(sum(abs(maha(X, mu, mcov) - bench)), 1e-6)
  expect_less_than(sum(abs(maha(X, mu, myChol, isChol = TRUE) - bench)), 1e-6)
  # Parallel
  expect_less_than(sum(abs(maha(X, mu, mcov, ncores = 2) - bench)), 1e-6)
  expect_less_than(sum(abs(maha(X, mu, myChol, ncores = 2, isChol = TRUE) - bench)), 1e-6)
  
  #############
  ###### d = 10, n = 10 case
  #############
  N <- 100
  d <- 10
  mu <- 1:d
  X <- t(t(matrix(rnorm(N*d), N, d)) + mu)
  tmp <- matrix(rnorm(d^2), d, d)
  mcov <- tcrossprod(tmp, tmp)
  myChol <- chol(mcov)
  
  ##### dmvn()
  bench <- dmvnorm(X, mu, mcov, log = T)
  # Sequential
  expect_less_than(sum(abs(dmvn(X, mu, mcov, log = T) - bench)), 1e-10)
  expect_less_than(sum(abs(dmvn(X, mu, myChol, isChol = TRUE, log = T) - bench)), 1e-10)
  # Parallel
  expect_less_than(sum(abs(dmvn(X, mu, mcov, ncores = 2, log = T) - bench)), 1e-10)
  expect_less_than(sum(abs(dmvn(X, mu, myChol, ncores = 2, isChol = TRUE, log = T) - bench)), 1e-10)
  
  ##### maha()
  bench <- mahalanobis(X, mu, mcov)
  # Sequential
  expect_less_than(sum(abs(maha(X, mu, mcov) - bench)), 1e-6)
  expect_less_than(sum(abs(maha(X, mu, myChol, isChol = TRUE) - bench)), 1e-6)
  # Parallel
  expect_less_than(sum(abs(maha(X, mu, mcov, ncores = 2) - bench)), 1e-6)
  expect_less_than(sum(abs(maha(X, mu, myChol, ncores = 2, isChol = TRUE) - bench)), 1e-6)
})