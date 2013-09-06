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

test_that("Test for monotonic regression function", {
  ## TODO: add more tests?
  A <- matrix(c(1, 0, 0, 1, 0, 1, 0, 0, 1), nrow=3)
  y <- c(2, 4, 8)

  output <- MR(y=y, E=A)

  x <- c(2, 6, 6)
  fit <- 8
  
  expect_that(x, equals(output$x))
  expect_that(fit, equals(output$fit))
})
