library(testthat)

source('../STA.R')

context("Conversion function tests")

test_that("Test for converting adjacency matrix to inequality coefficient matrix", {
  A <- matrix(c(1, 0, 0, 1, 0, 1, 0, 0, 1), nrow=3)
  adj2ineqA <- matrix(c(-1, 0, 0, 1, -1, 0, 0, -1, 1, 0, 0, -1), nrow=4, byrow=TRUE)

  expect_that(Adj2Ineq(A), equals(adj2ineqA))
})

test_that("Test for converting adjacency matrix to cell array", {
  A <- matrix(c(1, 0, 0, 1, 0, 1, 0, 0, 1), nrow=3)
  adj2cellA <- list(c(1, 1), c(1, 2), c(3, 2), c(3, 3))

  expect_that(Adj2Cell(A), equals(adj2cellA))
})

test_that("Test for converting cell array to adjacency matrix", {
  A <- list(c(1, 1), c(1, 2), c(3, 2), c(3, 3))
  cell2adjA <-  matrix(c(1, 0, 0, 1, 0, 1, 0, 0, 1), nrow=3)

  expect_that(cell2adjA, equals(Cell2Adj(1:3,A)))
})

context("Monotonic Regression tests")

test_that("Test for monotonic regression function with made up data", {
  A <- matrix(c(1, 0, 0, 1, 0, 1, 0, 0, 1), nrow=3)
  y <- c(2, 4, 8)

  output <- MR(y=y, E=A)

  x <- c(2, 6, 6)
  fit <- 8
  
  expect_that(x, equals(output$x))
  expect_that(fit, equals(output$fit))
})

## TODO: add more tests for MR with John's datasets?

context("State trace analysis statistics tests")

test_that("Test for state trace analysis statistics with x.dat", {
  x.dat <- matrix(scan('../../data/x.dat'), ncol = 5, byrow = TRUE);

  output <- staSTATS(x.dat)
  output <- output[[1]]
  
  means <- c(-0.06410421, 0.54172488, 1.25611257, 1.58557776, 2.13484796)
  n <- 20
  cov <- matrix(c(0.85214702, -0.08505141,  0.02161435,  0.06943527,  0.02321054,
                  -0.08505141,  1.02471854,  0.16462394, -0.02683712, -0.38845463,
                  0.02161435,  0.16462394,  1.27805730,  0.11634497, -0.08167361,
                  0.06943527, -0.02683712,  0.11634497,  1.72060522,  0.61117765,
                  0.02321054, -0.38845463, -0.08167361,  0.61117765,  2.93684274), nrow=5, byrow=TRUE)
  
  weights <- matrix(c(23.47013, 19.51755, 15.64875, 11.62382, 6.810034), nrow=1)
  
  expect_that(means, equals(output$means))
  expect_that(n, equals(output$n))
  expect_that(cov, equals(output$cov))
  expect_equal(weights, output$weights, tolerance=0.0001)
})

context("State trace analysis monotonic regression tests")

test_that("Test for state trace analysis monotonic regression", {
  x.dat <- matrix(scan('../../data/x.dat'), ncol = 5, byrow = TRUE);

  output <- staMR(x.dat)

  x <- c(-0.06410421, 0.54172488, 1.25611257, 1.58557776, 2.13484796)
  f <- 0

  expect_that(x, equals(output$x))
  expect_that(f, equals(output$f))
})
