library(pracma)
library(expm)

MR1 <- function(y, w, E) {
  # Uses lsqlin to solve standard monotonic regression problems
  # MR conducts monotonic regression on each column separately
  # y is a vector of values - the data
  # w is the weight matrix for the (lsqlin) fit: either a vector of positive
  # numbers or a positive definite matrix
  # E is an adjacency matrix coding the partial order model

  n <- length(y)

  if(nargs() == 1) {
    w <- diag(n)
    E <- matrix(0, n, n)
  }
  else if(nargs() == 2) {
    E <- matrix(0, n, n)
  }

  ## print(w)
  ## if(nrow(w) == 0 & ncol(w) == 0) {
  ##   w <- diag(n)
  ## }

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
    C <<- diag(sqrt(w))
  }
  else {
    C <<- expm::sqrtm(w) # use expm version of this function
  }

  d <<- C %*% y

  # TODO; find out what is wrong!
  
  (L <- lsqlin(A, b, C, d))
  return(C %*% L)
}

adj2cell <- function(adj) {
  # converts adjacency matrix to cell array
  
  row <- which(A != 0, arr.ind = TRUE)[, 1] # get row indices of non-zero elements
  column <- which(A != 0, arr.ind = TRUE)[, 2] # get column indices of non-zero elements

  E <- rep(list(list()), length(row)) # initializes a list of empty lists

  # fill in each list with corresponding row and column indices
  for(i in 1:length(row)) {
    E[[i]] <- c(row[i], column[i])
  }
  return(E)
}

adj2ineq <- function(adj) {
  # converts an adjacency matrix to an inequality coefficient matrix for monotonic regression
  
  i <- which(A != 0, arr.ind = TRUE)[, 1] # get row indices of non-zero elements
  j <- which(A != 0, arr.ind = TRUE)[, 2] # get column indices of non-zero elements
  s <- t(A)[t(A) != 0] # get values of non-zero elements, need to transpose to get correct order

  # initialize inequality matrix with zeros
  Aineq <- matrix(0, nrow = length(i), ncol = nrow(adj)) 

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

# for testing
A <- matrix(c(1, 1, 0, 0, 1, 1, 1, 1, 1), nrow = 3, byrow = TRUE)

E <- adj2cell(A)

nodes <- 1:nrow(A)

C <- cell2adj(nodes, E)

L <- MR1(c(1, 4, 3), c(1, 1, 1), A)
