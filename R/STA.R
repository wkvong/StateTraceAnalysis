library(limSolve)
library(expm)
library(R.matlab)
library(MASS)
library(pracma)

staPLOT <- function(data = c(), model = c(), groups = c(), labels = c(), axislabels = c(), axislimits = c()) {
  if(is.list(data)) {
    ys <- staSTATS(data)
  }
  else {
    ys <- outSTATS(data)
  }

  print(ys)
  
  x <- ys[[1]]$means
  y <- ys[[2]]$means

  ## TODO: Get different error bars depending on structure of ys
  cx <- sqrt(diag(ys[[1]]$cov)/diag(ys[[1]]$lm))
  cy <- sqrt(diag(ys[[2]]$cov)/diag(ys[[2]]$lm))

  if(isempty(groups)) {
    groups <- 1:length(ys[[1]]$means)
  }

  if(isempty(labels)) {
    labels <- c()
    for(i in 1:length(groups)) {
      labels[i] <- paste0("Condition ", i)
    }
  }

  if(isempty(axislabels)) {
    axislabels <- c()
    for(i in 1:length(groups)) {
      axislabels[i] <- paste0("Condition ", i)
    }
  }

  ## TODO: finish off! need to add different group conditions, titles/axes/error bars/axis limits
  plot(x,y)
}

CMRfits <- function(nsample, data, E = list(), E1 = list()) {
  ## TODO: function documentation
  
  if (is.list(data)) {
    type <- 0
    nvar <- length(unique(data[, 3]))
  }
  else {
    type <- 1
    nvar <- length(data)
  }

  if (type == 0) {
    ys <- staSTATS(data)
  }
  else {
    ys <- outSTATS(data)
  }

  if (isempty(E1) ) {
    if (!isempty(E)) {
      staMR.output <- staMR(ys, E)
      x2 <- staMR.output$x
      f2 <- staMR.output$f
    }
    else {
      f2 <- 0
    }

    staCMR.output <- staCMR(ys, E)
    x1 <- staCMR.output$x
    f1 <- staCMR.output$f
  }
  else {
      staCMR.output <- staCMR(ys, E)
      x2 <- staCMR.output$x
      f2 <- staCMR.output$f
      staCMR.output <- staCMR(ys, E1)
      x1 <- staCMR.output$x
      f1 <- staCMR.output$f
    }

  f <- f1 - f2
  datafit <- c(f, sum(f))

  ## set.generator("MersenneTwister", initialization="init2002", resolution=32, seed=12345)
  fits <- matrix(0, nsample, nvar)
  
  ## initiate parallel code
  ## cl <- makeCluster(4)
  ## registerDoParallel(cl)

  for (i in 1:nsample) {
    ## bootstrap sample
    yb <- bootstrap(data, type)
    ## fit 1D model to bootstrap data
    if (type == 0) {
      y <- staSTATS(yb)
    }
    else {
      y <- outSTATS(yb)
    }

    if (isempty(E1)) {
      staCMR.output <- staCMR(y, E)
      x <- staCMR.output$x
      f <- staCMR.output$f
    }
    else {
      staCMR.output <- staCMR(y, E1)
      x <- staCMR.output$x
      f <- staCMR.output$f
    }
    ## resample model
    yr <- resample(x, y, type)
    y <- yr

    if (isempty(E1)) {
      if (!isempty(E)) {
        staMR.output <- staMR(y, E)
        x2 <- staMR.output$x
        f2 <- staMR.output$f
      }
      else {
        f2 <- 0
      }

      staCMR.output <- staCMR(y, E)
      x1 <- staCMR.output$x
      f1 <- staCMR.output$f
    }
    else {
      staCMR.output <- staCMR(y, E)
      x2 <- staCMR.output$x
      f2 <- staCMR.output$f
      staCMR.output <- staCMR(y, E1)
      x1 <- staCMR.output$x
      f1 <- staCMR.output$f
    }
    f <- f1 - f2
    fits[i, ] <-  f ## store Monte Carlo fits
  }

  ## stop cluster
  ## stopCluster(cl)

  fits <- cbind(fits, rowSums(fits))
  
  p <- rep(0, length(datafit)) ## calculate p

  for (i in 1:ncol(fits)) {
    k <- which(fits[, i] >= datafit[i])
    p[i] <- length(k)/nsample
  }

  output <- list(p=p, datafit=datafit, fits=fits)
  return(output)
}

bootstrap <- function(y, type) {
  ## Draws bootstrap from data

  yb <- list()

  if (type == 0) {
    ## y in specific nsub x ncond format

    for (ivar in 1:length(y)) {
      a <- y[[ivar]]
      nsub <- nrow(a)
      yy <- matrix(0, nrow(a), ncol(a))
      b <- rep(0, nsub)
      v <- nsub

      while (v == nsub) {
        for (isub in 1:nsub) {
          u <- sample(nsub)
          yy[isub, ] <- a[u[1], ]
          b[isub] <- u[1]
        }

        v <- sum(b == b[1])
      }
      yb[[ivar]] <- yy
    }
  }
  else {
    ## y is in the general format
    
    cond <- sort(unique(y[, 2]))
    var <- sort(unique(y[, 3]))

    yb <- y

    for(j in 1:length(var)) {
      for(i in 1:length(cond)) {
        k <- which(y[, 2] == cond[i] & y[, 3] == var[j])
        a <- y[k, ]
        r <- floor(runif(length(k))*length(k))+1
        yb[k, ] <- a[r, ]
      }
    }
  }

  return(yb)
}

resample <- function(x, y, type) {
  yr <- y

  for (ivar in 1:length(x)) {
    if (type == 0) {
      sigma <- y[[ivar]]$cov/y[[ivar]]$n
    }
    else {
      sigma <- matrix(0, nrow(y[[ivar]]$cov), ncol(y[[ivar]]$cov))
      k <- which(y[[ivar]]$n > 0)
      sigma[k] <- y[[ivar]]$cov[k]/y[[ivar]]$n[k]
    }

    yr[[ivar]]$means <- mvrnorm(mu=x[[ivar]], Sigma=sigma)
  }

  return(yr)
}

staCMR <- function(data, E = list()) {
  ## Performs compound monotonic regression for state trace analysis
  ##
  ## Args:
  ##   data: List of submatrices containing data
  ##   E: List containing the starting partial order model
  ##
  ## Returns:
  ##   TODO
  
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

  CMR.output <- CMR(y, E)
  x <- CMR.output$x.star
  f <- CMR.output$f.fits
  e.prime <- CMR.output$e.star
  
  f[f < tol] <- 0

  if (!is.list(data)) {
    x <- x[[1]]
    f <- f[1]
  }

  output <- list(x=x, f=f, e.prime=e.prime)
  return(output)
}

CMR <- function(y, E = list()) {
  ## Coupled/Compound Monotonic Regression
  ## 
  ## Args:
  ##   y: A list of structured output from staSTATS (ie means, weights etc) 
  ##   E: A list of the starting partial order model: E for Edges
  ##      for example, if we want x3 <= x2 <= x1 and x4 <= x3 then
  ##      E = list(c(3, 2, 1), c(4, 3)) and so on.
  ##
  ## Returns a list with the following items:
  ##   x.star: Best fitting values
  ##   f.bar: Weighted least squares fit
  ##   e.star: Final partial order model 
  ##   exitflag: Vector of exit flags for fits for each variable

  L <- list()

  if (is.list(E)) {
    L[[1]] <- list()
    L[[1]][[1]] <- Cell2Adj(1:length(y[[1]]$means), E) ## E
  }
  else {
    L[[1]] <- list()
    L[[1]][[1]] <- E
  }

  L[[1]][[2]] <- -Inf ## F

  f.bar <- Inf
  e.bar <- L[[1]][[1]]

  exitflag <- list()
  
  while(length(L) > 0) {
    e.prime <- L[[1]][[1]]
    f.floor <- L[[1]][[2]]

    L[[1]] <- NULL ## Removes the first element from the list

    if (f.floor < f.bar) {
      staMR.output <- staMR(y, e.prime)
      x.prime <- staMR.output$x
      fits <- staMR.output$f
      exitflag <- staMR.output$exitflag
      
      f.fit <- sum(fits)

      if (f.fit < f.bar) {
        feas <- Feasible(x.prime) ## Check for feasible solution
        flag <- feas$flag
        idx <- feas$idx

        if (flag) {
          f.bar <- f.fit
          f.fits <- fits
          x.bar <- x.prime
          e.bar <- e.prime
        }
        else {
          ## Create two branches
          end <- length(L) + 1
          L[[end]] <- list()
          L[[end + 1]] <- list()

          L[[end]][[1]] <- e.prime
          L[[end]][[1]][idx[1], idx[2]] <- 1 ## (i,j) branch
          L[[end]][[2]] <- f.fit
          L[[end+1]][[1]] <- e.prime
          L[[end+1]][[1]][idx[2], idx[1]] <- 1 ## (j,i) branch
          L[[end+1]][[2]] <- f.fit
        }
      }
    }
  }

  x.star <- x.bar
  e.star <- e.bar

  output <- list(x.star=x.star, f.fits=f.fits, e.star=e.star, exitflag=exitflag)
  return(output)
}

Feasible <- function(x.prime) {
  ## Determines if there is a feasible solution (to what?)
  ##
  ## Args:
  ##   x.prime: List of possible solutions
  ##
  ## Returns:
  ##   True if there is a feasible solution, otherwise the largest inversion 
  
  flag <- 1
  idx <- vector()
  tol <- 1e-10
  x <- matrix(0, length(x.prime[[1]]), length(x.prime))
  
  for (i in 1:length(x.prime)) {
    x[, i] <- x.prime[[i]]
  }

  u <- combn(1:nrow(x), 2)
  u <- t(u) ## Transpose u to be consistent with matlab
  d <- x[u[, 1], ] - x[u[, 2], ]

  k <- which(abs(d) <= tol) 
  d[k] <- 0 ## Set zero to any tiny differences
  d <- sign(d) ## Convert differences to signs

  if (is.vector(d)) {
    s <- abs(d) - abs(d)
  }
  else {
    s <- rowSums(abs(d)) - abs(rowSums(d)) ## If signs of difference are not equal then s > 0
  }
    
  k <- which(s > 0)
  if (length(k) != 0) {
    flag <- 0
    idx <- u[k[1], ]
  }

  output <- list(flag=flag, idx=idx)
  return(output)
}

staMR <- function(data, E = list()) {
  ## Fits a monotonic regression mode for state trace analysis according
  ## to a specified partial order
  ##
  ## Args:
  ##   data: A list of data or structured output from staSTATS
  ##   E: A list or vector containing the partial ordering
  ##
  ## Returns a list with the following items:
  ##   x: Best fitting monotonic regression values to y-means
  ##   f: Fit statistic

  tol <- 10e-5

  if (!is.list(E) & is.vector(E)) {
    E <- list(E) ## Convert to list if vector is specified
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
  exitflag <- list()

  for(i in 1:length(y)) {
    y.i <- y[[i]]
    MR.output <- MR(y.i$means, y.i$weights, E)
    x[[i]] <- MR.output$x
    f[i] <- MR.output$fit
    exitflag[i] <- MR.output$exitflag
  }

  ## Round numbers close to tolerance level
  f[f < tol] <- 0

  if (!is.list(data)) {
    x <- x[[1]] ## TODO: check what this does?
  }

  output <- list(x=x, f=f, exitflag=exitflag)
  return(output)
}

outSTATS <- function(data) {
  ## Calculates statistics for state trace analysis for data in a "general" format
  
  ## Args:
  ##   data: A list of submatrices contain nsub x ncond matrices
  
  ## Returns a list with the following items:
  ##   means: Observed means
  ##   n: Number of subjects
  ##   cov: Observed covariance matrix
  ##   weights: Weight matrix for monotonic regression
  ##   lm: Loftus Masson statistics

  cond <- data[, 2]
  u.cond <- sort(unique(cond))
  var <- data[, 3]
  u.var <- sort(unique(var))
  within <- data[, 4:ncol(data)];

  ys <- list()
  
  for (i.var in 1:length(u.var)) {
    ys[[i.var]] <- list(c(), c(), c(), c(), c())
    names(ys[[i.var]]) <- c("means", "cov", "lm", "n", "weights")

    a <- list()
    b <- list()
    c <- list()
    
    for (i.cond in 1:length(u.cond)) {
      k <- which(var == u.var[i.var] & cond == u.cond[i.cond]);
      x <- within[k, ];
      u <- staSTATS(x);
      u <- u[[1]] ## staSTATS hack, needs to return everything in the first list element
      ys[[i.var]]$means <- c(ys[[i.var]]$means, u$means)
      a[[i.cond]] <- u$cov
      b[[i.cond]] <- u$lm
      c[[i.cond]] <- repmat(matrix(u$n), nrow(u$cov), ncol(u$cov))
    }

    s <- "ys[[i.var]]$cov <- blkdiag("

    for (i in 1:length(u.cond)) {
      s <- paste0(s, "a[[" , as.character(i), "]]")
      if (i < length(u.cond)) {
        s <- paste0(s, ", ")
      }
      else {
        s <- paste0(s, ")")
      }
    }

    eval(parse(text=s))

    s <- "ys[[i.var]]$lm <- blkdiag("

    for (i in 1:length(u.cond)) {
      s <- paste0(s, "b[[" , as.character(i), "]]")
      if (i < length(u.cond)) {
        s <- paste0(s, ", ")
      }
      else {
        s <- paste0(s, ")")
      }
    }

    eval(parse(text=s))

    s <- "ys[[i.var]]$n <- blkdiag("

    for (i in 1:length(u.cond)) {
      s <- paste0(s, "c[[" , as.character(i), "]]")
      if (i < length(u.cond)) {
        s <- paste0(s, ", ")
      }
      else {
        s <- paste0(s, ")")
      }
    }

    eval(parse(text=s))

    a <- diag(ys[[i.var]]$n)/diag(ys[[i.var]]$cov)
    ys[[i.var]]$weights <- t(a)
  }

  return(ys)
}

staSTATS <- function(data) {
  ## Calculates statistics for state trace analysis
  ##
  ## Args:
  ##   data: A list of submatrices contain nsub x ncond matrices
  ##
  ## Returns a list with the following items:
  ##   means: Observed means
  ##   n: Number of subjects
  ##   cov: Observed covariance matrix
  ##   weights: Weight matrix for monotonic regression
  ##   lm: TODO
  
  y <- data

  if (!is.list(data)) {
    y <- list(y)
  }

  output <- rep(list(list()), length(y))

  for(i in 1:length(y)) {
    y.i <- y[[i]]
    ## y.i <- y.i[complete.cases(y.i), ] ## delete rows with NAs

    out <- list()
    out.names <- c("means", "n", "cov", "weights", "lm")

    out[[1]] <- colMeans(y.i)
    out[[2]] <- nrow(y.i)

    if (out[[2]] > 1) {
      out[[3]] <- cov(y.i)
      out[[4]] <- nrow(y.i)/diag(cov(y.i))
    }
    else {
      out[[3]] <- matrix(0, nrow=ncol(y.i), ncol=ncol(y.i))
      out[[4]] <- matrix(1, nrow=ncol(y.i), ncol=1)
    }

    ## Tranpose weight matrix
    out[[4]] <- t(out[[4]])

    ## Calculate Loftus Masson data
    out[[5]] <- LoftusMasson(y.i)
    
    ## Add variable names to list
    names(out) <- out.names

    ## Add to output
    output[[i]] <- out
  }

  return(output)
}

LoftusMasson <- function(y.i) {
  y.cond <- repmat(matrix(colMeans(y.i), ncol=ncol(y.i)), nrow(y.i), 1)
  y.subj <- repmat(matrix(rowMeans(y.i), nrow=nrow(y.i)), 1, ncol(y.i))
  y.mean <- repmat(matrix(mean(rowMeans(y.i)),nrow=1), nrow(y.i), ncol(y.i))
  y.a <<- y.i - y.cond - y.subj + y.mean
  ss <- sum(y.a*y.a)
  df <- (nrow(y.i)-1)*(ncol(y.i)-1)
  ms.resid <- ss/df
  r <- diag(ms.resid, ncol(y.i))

  return(r)
}

repmat <- function(X, m, n) {
  ## R equivalent of repmat (matlab)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X, mx, nx*n)), mx*m, nx*n, byrow=TRUE)
}

MR <- function(y, w = diag(length(y)), E = matrix(0, length(y), length(y))) {
  ## Uses lsqlin to solve a standard monotonic regression problem
  ## MR conducts monotonic regression on each column separately
  ##
  ## Args:
  ##   y: Vector of values - the data
  ##   w: Weight matrix for the (lsqlin) fit: either a vector of positive
  ##      numbers or a positive definite matrix
  ##   E: Adjacency matrix coding the partial order model
  ##
  ## Returns a list containing the following:
  ##   x: Best fitting values
  ##   resid: TODO
  ##   fit: Weighted least squares fit
  ##   exitflag: Vector of exit flags for fits for each variable
  ##   output: Type of optimization used

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
    adj <- Cell2Adj(1:n, E)
  }
  else {
    adj <- E
  }

  if (sum(adj) > 0) {
    A <- -Adj2Ineq(adj) ## Turn adjacency matrix into a set of inequalities
    b <- matrix(0, nrow(A), 1)
  }
  else {
    A <- NULL ## not needed
    b <- NULL
  }

  if (is.vector(w) | nrow(w) == 1) {
    w <- w[1, ]
    C <- diag(sqrt(w))
  }
  else {
    C <- expm::sqrtm(w) ## Use expm version of this function
  }

  d <- C %*% y

  ## Run weighted least squares regression with inequality constraints
  L <- lsei(A = C, B = d, G = A, H = b, type=2, verbose = TRUE) 
  names(L) <- c("x", "resid", "fit", "exitflag", "output")

  return(L)
}

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
