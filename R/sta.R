setwd('/home/wkvong/Dropbox/research/current-projects/StateTraceAnalysis/R/')
source('stats.R')
source('convert.R')

library(limSolve)
library(expm)
library(R.matlab) ## TODO: remove eventually
library(MASS)
library(pracma)
library(reshape2)
library(plyr)
library(ggplot2)

CMRMNfits <- function(nsample, data, E = list()) {
  ## Parametric bootstrap sampling test of multinomial data
  ## data is NCOND x NRESP matrix of counts
  ## E is partial order in form of cell array or adjacency matrix

  y <- data
  staMRMN.output <- staMRMN(y, E)
  x <- staMRMN.output$x
  f2 <- staMRMN.output$f2
  staCMRMN.output <- staCMRMN(y, E)
  x <- staCMRMN.output$x
  f1 <- staCMRMN.output$f1

  datafit <- c(f1, f2, f1-f2)

  fits <- matrix(0, nsample, 3)

  for (i in 1:nsample) {
    ## bootstrap sample
    yb <- resampleMN(y)

    ## fit model to bootstrap sample
    staCMRMN.output <- staCMRMN(yb, E)
    x <- staCMRMN.output$x

    ## resample model
    yr <- resampleMN(x)

    ## fit model
    staMRMN.output <- staMRMN(yr, E)
    x <- staMRMN.output$x
    f2 <- staMRMN.output$fits ## TODO: check if f or fits is returned

    staCMRMN.output <- staCMRMN(yr, E)
    x <- staCMRMN.output$x
    f1 <- staCMRMN.output$f

    fits[isample: ] <- c(f1, f2, f1-f2)
  }

  ## calculate p
  k <- which(fits[, 3] >= datafit[3])
  p <- length(k)/nsample

  output <- list(p=p, datafit=datafit, fits=fits)
  return(output)
}

resampleMN <- function(x) {
  s <- rowSums(x)
  n <- round(s)
  p <- x/repmat(s, 1, ncol(x))
  yr <- rmultinom(1, n, p)

  return(yr)
}

MRMNfits <- function(nsample, data, E) {
  
}

staCMRMN <- function(y, E = list(), flag = 0) {
  L <- list()

  if (is.list(E)) {
    L[[1]][[1]] <- Cell2Adj(1:nrow(y[[1]], E)) ## E
  }
  else {
    L[[1]][[1]] <- E
  }

  L[[1]][[2]] <- -Inf ## F

  f.bar <- Inf
  e.bar <- L[[1]][[1]]

  while(length(L) > 0) {
    e.prime <- L[[1]][[1]]
    f.floor <- L[[1]][[2]]

    L[[1]] <- NULL ## Removes the first element from the list

    if (f.floor < f.bar) {
      staMRMN1.output <- staMRMN1(y, e.prime, flag) ## Least squares optimization
      x.prime <- staMRMN1.output$x
      fit <- staMRMN1.output$fit
      exitflag$ staMRMN1.output$exitflag

      f.fit <- fit

      if (f.fit < f.bar) {
        feas <- FeasibleCMRMN(x.prime) # Check if feasible, else get an index
        flag <- feas$flag
        idx <- feas$idx

        if (flag) {
          f.bar <- f.fit
          x.bar <- x.prime
          e.bar <- e.prime
        }
        else {
          end <- length(L) + 1
          L[[end]] <- list()
          L[[end+1]]
          
          L[[end]][[1]] <- e.prime
          L[[end]][idx[1], idx[2]] <- 1 ## (i,j) branch
          L[[end]][[2]] <- f.fit
          L[[end+1]][[1]] <- e.prime
          L[[end+1]][idx[2], idx[1]] <- 1 ## (j, i) branch
          L[[end+1]][[2]] <- f.fit
        }
      }
    }
  }

  x.star <- x.bar
  e.star <- e.bar

  k <- which(x.star < 0)
  x.star[k] <- 0 ## round negative values to 0
  g2.fit <- NaN

  output <- list(x.star=x.star, f.bar=f.bar, g2.fit=g2.fit, e.star=e.star)
  return(output)
}

testCMRMN <- function() {
  y <- matrix(c(66, 16, 2, 3, 7, 11, 102, 14, 6, 7, 4, 5, 94, 12, 7, 4, 4, 6, 34, 10, 21, 13, 18, 22, 46, 19, 20, 12, 8, 19, 58, 22, 12, 13, 13, 13, 36, 8, 10, 14, 15, 24, 47, 11, 10, 11, 13, 34, 48, 7, 8, 8, 13, 49, 12, 9, 13, 22, 19, 44, 24, 16, 13, 24, 21, 44, 19, 14, 18, 20, 26, 37), nrow=12, byrow=TRUE)
  E <- list(c(1, 2, 3), c(4, 5, 6), c(4, 1), c(5, 2), c(6, 3), c(7, 8, 9), c(10, 11, 12), c(10, 7), c(11, 8), c(12, 9))
  staCMRMN(y, E)
}

FeasibleCMRMN <- function(x.prime) {
  flag <- 1
  idx <- c()
  tol <- 1e-10
  
  
  output <- list(flag=flag, idx=idx)
  return(output)
}

staMRMN1 <- function(y, e.prime, flag) {
  n.sum <- repmat(ncol(y), 1, ncol(y))
  p <- p/n.sum
  dim(p) <- c(length(p), 1)
  dim(n.sum) <- c(length(n.sum, 1))
  weights <- diag(as.vector(n.sum))

  ## augment for all columns
  a <- Adj2Ineq(E)
  A <- ineqrep(a, ncol(y))
  b <- matrix(0, nrow(A), 1)

  Aeq <- matrix()
  for(j = 1:ncol(y)) {
    Aeq <- cbind(Aeq, matrix(1, nrow(y), nrow(y)))
  }

  beq <- matrix(1, nrow(y), nrow(y))

  ## Set up matrices and bounds
  C <- expm::sqrtm(weights)
  d <- C %*% p
  x0 <- p
  lb <- matrix(0, nrow(x0), ncol(x0))
  ub <- matrix(1, nrow(x0), ncol(x0))

  L <- lsei(C, d, A, b, Aeq, beq) ## TODO: figure out how to do this properly!
  x1 <- L$x

  ## finish off with ML optimization
  if(L$flag > 0) {
    ## TODO: find R equivalent of fmincon here!!
  }
}

ML <- function(x, y) {
  if (is.vector(x)) {
    ## if x is a vector of probabilities
    dim(y) <- c(length(y), 1)
    k <- which(x == 0)
    p <- x
    p[k] <- 1/sum(y)
    u <- y*log(p)
    m <- -sum(u)
  }
  else {
    ## if x is a matrix of expected counts
    p <- x/repmat(ncol(y), 1, ncol(y))
    k <- which(x == 0)
    p[k] <- 1/sum(y)
    u <- y*log(p)
    m <- -sum(sum(u))
  }
}

ineqrep <- function(a, n) {
  b <- matrix()

  for (i in 1:n) {
    u <- c()
    for (j in 1:n) {
      if (j <= i) {
        u <- cbind(u, a)
      }
      else {
        u <- cbind(u, matrix(0, length(a), length(a)))
      }
    }
    b <- rbind(b, u)
  }

  return(b)
}

staMRMN <- function(y, E, flag) {
  
}

CMRfits <- function(nsample, data, E = list(), E1 = list()) {
  ## TODO: function documentation

  ## TODO: ask John about determining correct structure type here
  if (is.list(data)) {
    type <- 0
    if (length(data) > 2) {
      nvar <- length(unique(data[, 3]))
    }
    else {
      nvar <- length(data)
    }
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

    for (j in 1:length(var)) {
      for (i in 1:length(cond)) {
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
  
  if (!is.list(data)) {
    y <- outSTATS(y)
  }
  else if (is.list(data[[1]])) {
    y <- data
  }
  else {
    y <- staSTATS(y)
  }
  
  CMR.output <- CMR(y, E)
  x <- CMR.output$x.star
  f <- CMR.output$f.fits
  e.prime <- CMR.output$e.star
  
  f[f < tol] <- 0

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
        feas <- FeasibleCMR(x.prime) ## Check for feasible solution
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

FeasibleCMR <- function(x.prime) {
  ## Determines if there is a feasible solution (to what?)
  ##
  ## Args:
  ##   x.prime: List of possible solutions
  ##
  ## Returns:
  ##   True if there is a feasible solution, otherwise the largest inversion 
  
  flag <- 1
  idx <- c()
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
    s <- sum(abs(d))-abs(sum(d))  #this is a John Fix
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

  for (i in 1:length(y)) {
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
