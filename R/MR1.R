## TODO: give names to lists and call them when appropriate

library(limSolve)
library(expm)

staMR <- function(data, E = list()) {
  tol <- 10e-5

  if(!is.list(E) & is.vector(E)) {
    E <- list(E) ## convert to list if vector is specified
  }

  y <- data

  if(!is.list(y)) {
    y <- list(y)
  }

  if(!is.list(y[[1]])) {
    y = staSTATS(y);
  }

  x <- list()
  f <- rep(0, length(y))
  ## TODO: exitflag information

  names(y) <- c("means", "weights")
  
  for(i in 1:length(y)) {
    ## TODO: y[[i]].means/weights
    reg <- MR1(y[[i]], y[[i]], E)

    
    x[[i]] <- reg$X
    f[i] <- reg$solutionNorm
    ## TODO: exit flag
  }

  ## round numbers close to tolerance level
  f[f < tol] <- 0

  if(!is.list(data)) {
    x <- x[[1]] ## TODO: what does this do?
  }

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

  if(!is.list(data)) {
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

    if(out[[2]] > 1) {
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

    print(out)
    
    ## add to output
    output[[idata]] <- out
  }

  if(!is.list(data)) {
    output <- output[[1]]
  }

  return(output)
}

MR1 <- function(y, w = diag(length(y)), E = matrix(0, length(y), length(y))) {
  # Uses lsqlin to solve standard monotonic regression problems
  # MR conducts monotonic regression on each column separately
  # y is a vector of values - the data
  # w is the weight matrix for the (lsqlin) fit: either a vector of positive
  # numbers or a positive definite matrix
  # E is an adjacency matrix coding the partial order model

  n <- length(y)

  if(nrow(w) == 0 & ncol(w) == 0) {
    w <- diag(n)
  }

  if(nrow(E) == 0 & ncol(E) == 0) {
    E <- matrix(0, n, n)
  }

  if(is.list(E)) {
    adj <- cell2adj(1:n, E)
  }
  else {
    adj <- E
  }

  if(sum(adj) > 0) {
    A <- adj2ineq(adj) # turn adjacency matrix into a set of inequalities
    b <- matrix(0, nrow(A), 1)
  }
  else {
    A <- matrix()
    b <- matrix()
  }

  if(is.vector(w)) {
    C <- diag(sqrt(w))
  }
  else {
    C <- expm::sqrtm(w) # use expm version of this function
  }

  d <- C %*% y

  L <- lsei(A = C, B = d, G = -A, H = b, type=2, verbose = TRUE)
  return(L)
}

adj2cell <- function(adj) {
  # converts adjacency matrix to cell array
  
  row <- which(adj != 0, arr.ind = TRUE)[, 1] # get row indices of non-zero elements
  column <- which(adj != 0, arr.ind = TRUE)[, 2] # get column indices of non-zero elements

  E <- rep(list(list()), length(row)) # initializes a list of empty lists

  # fill in each list with corresponding row and column indices
  for(i in 1:length(row)) {
    E[[i]] <- c(row[i], column[i])
  }
  return(E)
}

adj2ineq <- function(adj) {
  # converts an adjacency matrix to an inequality coefficient matrix for monotonic regression
  
  i <- which(adj != 0, arr.ind = TRUE)[, 1] # get row indices of non-zero elements
  j <- which(adj != 0, arr.ind = TRUE)[, 2] # get column indices of non-zero elements
  s <- t(adj)[t(adj) != 0] # get values of non-zero elements, need to transpose to get correct order

  print(i)
  print(j)
  
  # initialize inequality matrix with zeros
  Aineq <- matrix(0, nrow = length(i), ncol = nrow(adj)) 

  print(Aineq)
  
  # fill in matrix with inequality coefficients
  for(k in 1:length(i)) {
    Aineq[k, i[k]] <- 1;
    Aineq[k, j[k]] <- -1;
  }

  return(Aineq)
}

cell2adj <- function(nodes, E=list()) {
  # converts a partial order model in cell array form to an adjacency matrix suitable for monotonic regression
  
  if(!is.list(E))
    E <- list(E) # probably won't work unless input is a set of lists, TODO: ask john about this

  n <- length(nodes)
  
  # initialize adjacency matrix with zeros
  adj <- matrix(0, nrow = n, ncol = n)

  # fill in adjacency matrix
  if(length(E) != 0) {
    for(i in 1:length(E)) {
      if(length(E[[i]]) != 0) {
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

A <- matrix(c(1, 0, 1, 0, 1, 0, 0, 1, 0), nrow = 3, byrow = TRUE)

w <- diag(3)

L <- MR1(y=c(2, 4, 8), w, E=A)

print(L)

## read in data and convert to a list
x <- matrix(scan('../matlab/x.dat'), ncol = 5, byrow = TRUE)

output <- staSTATS(x)
print(output)
