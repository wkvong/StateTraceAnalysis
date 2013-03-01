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

  Aineq <- matrix(0, length(i), nrow(adj)) # create zero matrix 

  # fill in matrix with inequality coefficients
  for(k in 1:length(i)) {
    Aineq[k, i[k]] = 1;
    Aineq[k, j[k]] = -1;
  }

  return(Aineq)
}

cell2adj <- function(nodes, E=list()) {
  return(E)
}

# for testing
A <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2, byrow = TRUE)

A <- matrix(c(1, 2, 0, 0, 2, 1, 1, 2, 3), nrow = 3, byrow = TRUE)

B <- adj2ineq(A)

C <- cell2adj(c(1, 2))
