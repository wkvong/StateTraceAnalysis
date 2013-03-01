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
    Aineq[k, i[k]] = 1;
    Aineq[k, j[k]] = -1;
  }

  return(Aineq)
}

cell2adj <- function(nodes, E=list()) {
  # converts a partial order model in cell array form to an adjacency matrix suitable for monotonic regression
  
  if(!is.list(E))
    E <- list(E) # probably won't work unless input is a set of lists, TODO: ask john about this

  n <- length(nodes)
  print(n)
  
  # initialize adjacency matrix with zeros
  adj <- matrix(0, nrow = n, ncol = n)

  # fill in adjacency matrix
  if(length(E) != 0) {
    for(i in 1:length(E)) {
      if(length(E[[i]]) != 0) {
        # TODO: ask john what is wrong
        ## u <- choose(E[[i]], 2)
        u <- matrix(E[[i]], nrow = 1) # convert u to matrix to use nrow
        for(j in 1:nrow(u)) {
          k1 <- which(nodes == u[j, 1])
          k2 <- which(nodes == u[j, 2])
          adj[k1, k2] = 1;
        }
      }
    }
  }
  
  return(adj)
}

# for testing
A <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2, byrow = TRUE)

A <- matrix(c(1, 1, 0, 0, 1, 1, 1, 1, 1), nrow = 3, byrow = TRUE)

E <- adj2cell(A)

nodes <- c(1, 2, 3)

C <- cell2adj(nodes, E)
