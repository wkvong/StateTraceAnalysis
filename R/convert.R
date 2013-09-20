## Conversion functions between different data structures

Adj2Cell <- function(adj) {
  ## Converts adjacency matrix to list
  ##
  ## Args:
  ##   adj: Adjacency matrix
  ##
  ## Returns:
  ##   List
  
  row <- which(adj != 0, arr.ind = TRUE)[, 1] ## Get row indices of non-zero elements
  column <- which(adj != 0, arr.ind = TRUE)[, 2] ## Get column indices of non-zero elements

  E <- rep(list(list()), length(row)) ## Initializes a list of empty lists

  ## Fill in each list with corresponding row and column indices
  for(i in 1:length(row)) {
    E[[i]] <- c(row[i], column[i])
  }

  return(E)
}

Adj2Ineq <- function(adj) {
  ## Converts an adjacency matrix to an inequality coefficient matrix 
  ##
  ## Args:
  ##   adj: Adjacency matrix
  ##
  ## Returns:
  ##   Matrix containing inequality coefficients
  
  i <- which(adj != 0, arr.ind = TRUE)[, 1] ## Get row indices of non-zero elements
  j <- which(adj != 0, arr.ind = TRUE)[, 2] ## Get column indices of non-zero elements
  s <- t(adj)[t(adj) != 0] ## Get values of non-zero elements, need to transpose to get correct order

  ## Initialize inequality matrix with zeros
  Aineq <- matrix(0, nrow = length(i), ncol = nrow(adj)) 
  
  ## Fill in matrix with inequality coefficients
  for(k in 1:length(i)) {
    Aineq[k, i[k]] <- 1;
    Aineq[k, j[k]] <- -1;
  }

  return(Aineq)
}

Cell2Adj <- function(nodes, E=list()) {
  ## Converts a partial order model in list form to an adjacency matrix suitable for monotonic regression
  ##
  ## Args:
  ##   nodes: TODO?
  ##   E: List containing the partial order model
  ##
  ## Returns:
  ##   Adjacency matrix for monotonic regression
  
  if (!is.list(E))
    E <- list(E)

  n <- length(nodes)
  
  ## Initialize adjacency matrix with zeros
  adj <- matrix(0, nrow = n, ncol = n)

  ## Fill in adjacency matrix
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
