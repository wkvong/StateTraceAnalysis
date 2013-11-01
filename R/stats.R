## Calculate statistics

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
