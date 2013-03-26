## TODO: give names to lists and call them when appropriate

library(limSolve)
library(expm)

## TODO: function documentation
staCMR <- function(data, E = list()) {
  tol <- 10e-5

  if (!is.list(E)) {
    E <- list(E)
  }
  
  y <- data

  if (!is.list(y)) {
    y <- list(y)
  }

  if (!is.list(y[[1]])) {
    y <- staSTATS(y)
  }

  output <- CMR(y, E)

  ## TODO: set x/f/Eprime based on values from output
  x <- 0
  f <- 0
  Eprime <- 0
  
  f[f < tol] <- 0

  ## TODO: check this work
  if (!is.list(data)) {
    x <- x[[1]]
    f <- f[1]
  }

  return(list(x, f, Eprime))
}

CMR <- function(y, E) {
  ## TODO: clean up function documentation to R standards
  ## coupled monotonic regression
  ## y is a cell array of structured output from staSTATS (ie means, weights etc) 
  ## E is a cell array of the starting partial order model: E for Edges
  ## for example, if we want x3 <= x2 <= x1 and x4 <= x3 then
  ## E = {[3 2 1] [4 3]}, and so on.
  ## Fbar is the fit of the starting model -- 'Inf' is a good start
  ## returns:
  ## xystar = best fitting values
  ## Fbar = weighted least squares fit
  ## Estar = final partial order model 
  ## exitflag = vector of exit flags for fits for each variable

  L <- list()

  if (is.list(E)) {
    L[[1]] <- list()
    L[[1]][[1]] <- cell2adj(1:length(y[[1]]$means), E)
    names(L[[1]]) <- c("E")
  }
  else {
    L[[1]] <- list()
    L[[1]][[1]] <- E
    names(L[[1]]) <- c("E")    
  }

  L[[1]][[2]] <- -Inf
  Fbar <- Inf
  EBar <- L[[1]][[1]] ## TODO: rewrite as L[[1]]$E?

  while(length(L) > 0) {
    Eprime <- L[[1]][[1]] ## TODO: rewrite with variables names?
    Ffloor <- L[[1]][[2]]

    L[[1]] <- NULL ## remove this element from the list
    if (Ffloor < Fbar) {
      output <- staMR(y, Eprime)
      names(output) <- c("xPrime", "fits") ## TODO: add exitflag information
      Ffit <- sum(output$fits)

      if (Ffit < Fbar) {
        feas <- Feasible(output$xPrime)
        names(feas) <- c("flag", "idx")

        if (feas$flag) {
          Fbar <- Ffit
          Ffits <- fits
          xBar <- xPrime
          EBar <- Eprime
        }
        else {
          ## create two branches
          end <- length(L)
          L[[end + 1]] <- list()
          L[[end + 1]][[1]] <- Eprime
          ## TODO: finish this off
        }
      }
    }
  }

  xStar <- xBar
  EStar <- EBar

  ## TODO: return exitflag
  return(list(xStar, Ffits, EStar))
}

## TODO: test this function
Feasible <- function(xx) {
  flag <- 1
  idx <- vector()
  zero <- 1e-10
  x <- matrix(0, length(xx[[1]]), length(xx))

  for (i in 1:length(xx)) {
    x[, i] <- xx[[i]]
  }

  u <- combn(1:nrow(x), 2)
  d <- x[u[, 1], ] - x[u[, 2], ]
}

## TODO: function documentation
staMR <- function(data, E = list()) {
  tol <- 10e-5

  if (!is.list(E) & is.vector(E)) {
    E <- list(E) ## convert to list if vector is specified
  }

  y <- data

  if (!is.list(y)) {
    y <- list(y)
  }

  if (!is.list(y[[1]])) {

    y = staSTATS(y);
  }
  
  x <- list()
  f <- rep(0, length(y))
  ## TODO: exitflag information

  for(i in 1:length(y)) {
    yy <- y[[i]]

    output <- MR1(yy$means, yy$weights, E)
    
    x[[i]] <- output$X
    f[i] <- output$solutionNorm
    ## TODO: exit flag
  }

  ## round numbers close to tolerance level
  f[f < tol] <- 0

  ## if (!is.list(data)) {
  ##   x <- x[[1]] ## TODO: what does this do?
  ## }

  ## TODO: return exitflag as well
  return(list(x, f))
}

staSTATS <- function(data) {
  ## TODO: reformat function documentation to fit R standards
  ## data is NSUB x NCOND sub-matrix or cell array of sub-matrices
  ## returns means, cov, nsub and weights
  ## output.means = observed means
  ## output.n = number of subjects
  ## output.cov = observed covariance matrix
  ## output.weights = weight matrix for monotonic regression
  
  y <- data

  if (!is.list(data)) {
    y <- list(y)
  }

  output <- rep(list(list()), length(y))

  for(idata in 1:length(y)) {
    yy <- y[[idata]]
    yy <- yy[complete.cases(yy), ] ## delete rows with NAs

    out <- list()
    out.names <- c("means", "n", "cov", "weights", "lm")

    out[[1]] <- colMeans(yy)
    out[[2]] <- nrow(yy)

    if (out[[2]] > 1) {
      out[[3]] <- cov(yy)
      out[[4]] <- nrow(yy)/diag(cov(yy))
    }
    else {
      out[[3]] <- matrix(0, nrow=ncol(yy), ncol=ncol(yy))
      out[[4]] <- matrix(1, nrow=ncol(yy), ncol=1)
    }

    ## tranpose weight matrix
    out[[4]] <- t(out[[4]])
    out[[5]] <- diag(diag(out[[3]]))

    ## give names to out list
    names(out) <- out.names

    ## add to output
    output[[idata]] <- out
  }

  if (!is.list(output)) {
    output <- output[[1]]
  }

  return(output)
}

MR1 <- function(y, w = diag(length(y)), E = matrix(0, length(y), length(y))) {
  ## Uses lsqlin to solve standard monotonic regression problems
  ## MR conducts monotonic regression on each column separately
  ## y is a vector of values - the data
  ## w is the weight matrix for the (lsqlin) fit: either a vector of positive
  ## numbers or a positive definite matrix
  ## E is an adjacency matrix coding the partial order model

  n <- length(y)

  if (nrow(w) == 0 & ncol(w) == 0) {
    w <- diag(n)
  }

  if (is.matrix(E)) {
    if (nrow(E) == 0 & ncol(E) == 0) {
      E <- matrix(0, n, n)
    }
  }

  if (is.list(E)) {
    adj <- cell2adj(1:n, E)
  }
  else {
    adj <- E
  }

  if (sum(adj) > 0) {
    A <- -adj2ineq(adj) ## turn adjacency matrix into a set of inequalities
    b <- matrix(0, nrow(A), 1)
  }
  else {
    A <- NULL
    b <- NULL
  }

  if (is.vector(w) | nrow(w) == 1) {
    w <- w[1, ]
    C <- diag(sqrt(w))
  }
  else {
    C <- expm::sqrtm(w) ## use expm version of this function
  }

  d <- C %*% y

  L <- lsei(A = C, B = d, G = A, H = b, type=2, verbose = TRUE)

  return(L)
}

adj2cell <- function(adj) {
  ## converts adjacency matrix to cell array
  
  row <- which(adj != 0, arr.ind = TRUE)[, 1] ## get row indices of non-zero elements
  column <- which(adj != 0, arr.ind = TRUE)[, 2] ## get column indices of non-zero elements

  E <- rep(list(list()), length(row)) ## initializes a list of empty lists

  ## fill in each list with corresponding row and column indices
  for(i in 1:length(row)) {
    E[[i]] <- c(row[i], column[i])
  }
  return(E)
}

adj2ineq <- function(adj) {
  ## converts an adjacency matrix to an inequality coefficient matrix for monotonic regression
  
  i <- which(adj != 0, arr.ind = TRUE)[, 1] ## get row indices of non-zero elements
  j <- which(adj != 0, arr.ind = TRUE)[, 2] ## get column indices of non-zero elements
  s <- t(adj)[t(adj) != 0] ## get values of non-zero elements, need to transpose to get correct order

  ## initialize inequality matrix with zeros
  Aineq <- matrix(0, nrow = length(i), ncol = nrow(adj)) 
  
  ## fill in matrix with inequality coefficients
  for(k in 1:length(i)) {
    Aineq[k, i[k]] <- 1;
    Aineq[k, j[k]] <- -1;
  }

  return(Aineq)
}

cell2adj <- function(nodes, E=list()) {
  ## converts a partial order model in cell array form to an adjacency matrix suitable for monotonic regression
  
  if (!is.list(E))
    E <- list(E)

  n <- length(nodes)
  
  ## initialize adjacency matrix with zeros
  adj <- matrix(0, nrow = n, ncol = n)

  ## fill in adjacency matrix
  if (length(E) != 0) {
    for(i in 1:length(E)) {
      if (length(E[[i]]) != 0) {
        u <- t(combn(E[[i]], 2))
        for(j in 1:nrow(u)) {
          k1 <- which(nodes == u[j, 1])
          k2 <- which(nodes == u[j, 2])
          adj[k1, k2] <- 1;
        }
      }
    }
  }
  return(adj)
}

## testing code

A <- matrix(c(1, 0, 1, 0, 1, 0, 0, 1, 0), nrow = 3, byrow = TRUE)

w <- diag(3)

L <- MR1(y=c(2, 4, 8), w, E=A)

## read in data and convert to a list
x <- matrix(scan('../matlab/x.dat'), ncol = 5, byrow = TRUE)

output <- staSTATS(x)
## print(output)

staMR(x)

CMR(output, list(c(3, 2, 1), c(4, 3)))
