source('sta.R')

staPLOT <- function(data = c(), model = c(), groups = c(), labels = c(), xlab = c(), ylab = c(), xlim = c(), ylim = c()) {
  ## Generates a state-trace plot

  ## TODO: arguments and output
  
  if (is.list(data)) {
    ys <- staSTATS(data) ## within-subjects structured data
  }
  else {
    ys <- outSTATS(data) ## general data format
  }
  ## TODO: case if data is already is already in stats form
  
  x <- ys[[1]]$means
  y <- ys[[2]]$means

  ## TODO: Get different error bars depending on structure of ys
  cx <- sqrt(diag(ys[[1]]$cov)/diag(ys[[1]]$lm))
  cy <- sqrt(diag(ys[[2]]$cov)/diag(ys[[2]]$lm))

  if (isempty(groups)) {
    groups <- 1:length(ys[[1]]$means)
  }

  if (isempty(labels)) {
    labels <- c()
    for (i in 1:length(groups)) {
      labels[i] <- paste0("Condition ", i)
    }
  }

  if (isempty(xlab)) {
    xlab <- "Outcome Variable 1"
  }

  if (isempty(ylab)) {
    ylab <- "Outcome Variable 2"
  }

  ## TODO: plot points
  plotdata(x, y, groups, 0)

  ## TODO: plot error bars
  
  ## TODO: plot model
  if (!isempty(model)) {
    xm <- model[[1]]
    ym <- model[[2]]
    ix <- tiesort(xm, ym)
    ## TODO: plot model here
  }

  ## TODO: axis limits/axis labels

  ## TODO: legend
}

plotdata <- function(x, y, groups, flag) {

  for (i in 1:length(groups)) {
    a <- x[groups[[i]]]
    b <- y[groups[[i]]]
    print(a)
    print(b)
    
    if (i == 1) {
      print("print group 1")
      plot(a, b, xlim=c(0.3, 0.8), ylim=c(0.2, 0.7), pch=16)
    }
    if (i == 2) {
      print("print group 2")
      points(a, b)
    }
    ## TODO: other points
  }
  
  if (flag) {
    ## DeleteLegendEntry(h)
  }
}

tiesort <- function(x, y) {
  ## Sorts y values in increasing order
  ## Tied x values will sort y values in increasing order
  bignumber <- 100
  x <- round(x*bignumber)/bignumber
  y <- round(y*bignumber)/bignumber
  tt <- 1:length(x)
  z <- c(x, y, t(tt))
  ## TODO: sortrows
  a <- c()
  ix <- c()
  ## ix <- a[:3]

  return(ix)
}

testPlot <- function() {
  y <- matrix(scan('../data/delay.dat'), ncol = 7, byrow = TRUE)
  output <- staCMR(delay)
  x <- output$x
  groups <- list(1:4, 5:8)
  staPLOT(data=y, model=x, groups=groups)
}

