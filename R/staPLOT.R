source('sta.R')

staPLOT <- function(data = c(), model = c(), groups = c(), lab = c(), xlab = c(), ylab = c(), xlim = c(0, 1), ylim = c(0, 1)) {
  ## Generates a state-trace plot

  ## TODO: arguments and output documentation

  if (is.list(data)) {
    ys <- staSTATS(data) ## within-subjects structured data
  }
  else {
    ys <- outSTATS(data) ## general data format
  }
  ## TODO: cases if data is already is already in stats form
  
  x <- ys[[1]]$means
  y <- ys[[2]]$means

  ## TODO: Get different error bars depending on structure of ys
  if (is.list(data)) {
    cx <- sqrt(diag(ys[[1]]$cov)/ys[[1]]$n)
    cy <- sqrt(diag(ys[[2]]$cov)/ys[[2]]$n)
  }
  else {
    cx <- sqrt(diag(ys[[1]]$cov)/diag(ys[[1]]$n))
    cy <- sqrt(diag(ys[[2]]$cov)/diag(ys[[2]]$n))
  }

  if (isempty(groups)) {
    groups <- list()
    for(i in 1:length(ys[[1]]$means)) {
      groups[[i]] <- i
    }
  }

  if (isempty(lab)) {
    lab <- c()
    for (i in 1:length(groups)) {
      lab[i] <- paste0("Condition ", i)
    }
  }

  if (isempty(xlab)) {
    xlab <- "Outcome Variable 1"
  }

  if (isempty(ylab)) {
    ylab <- "Outcome Variable 2"
  }

  attr(groups, "varname") <- "groups"
  sta.df <- melt(groups)
  sta.df <- cbind(sta.df, x)
  sta.df <- cbind(sta.df, y)
  sta.df <- cbind(sta.df, cx)
  sta.df <- cbind(sta.df, cy)

  if (!isempty(model)) {
    xm <- model[[1]]
    ym <- model[[2]]
    ix <- tiesort(xm, ym)
    sta.df <- cbind(sta.df, xm=xm[ix])
    sta.df <- cbind(sta.df, ym=ym[ix])
  }

  sta.df <- sta.df

  cols <- c() ## TODO: add different colours here
  
  ggplot(sta.df, aes(x=x, y=y, colour=factor(groups))) +
    geom_errorbarh(aes(xmin=x-cx, xmax=x+cx), width=0.01) +
    geom_errorbar(aes(ymin=y-cy, ymax=y+cy), width=0.01) +
    geom_line(aes(x=xm, y=ym), linetype="dotted") +
    geom_point(aes(fill=factor(groups)), size=5, shape=21) +
    xlab(xlab) +
    ylab(ylab) +
    xlim(xlim) +
    ylim(ylim) +
    scale_colour_manual(values=c("black", "black"),
                        labels=lab) +
    scale_fill_manual(values=c("black", "white"),
                      labels=lab) +
    theme_bw() +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
    theme(legend.title=element_blank()) +
    theme(legend.justification=c(1,0)) +
    theme(legend.position=c(1,0))
}

tiesort <- function(x, y) {
  ## Sorts y values in increasing order
  ## Tied x values will sort y values in increasing order
  bignumber <- 100
  x <- round(x*bignumber)/bignumber
  y <- round(y*bignumber)/bignumber
  tt <- 1:length(x)
  z <- data.frame(cbind(x, y, tt))
  a <<- arrange(z, x, y)
  ix <- a[, 3]

  return(ix)
}

testPlot <- function() {
  y <- matrix(scan('../data/delay.dat'), ncol = 7, byrow = TRUE)
  output <- staCMR(y)
  x <- output$x
  groups <- list(1:4, 5:8)
  staPLOT(data=y, model=x, groups=groups, lab=c("No Delay", "Delay"), xlab="RB", ylab="II", xlim=c(0.3, 0.8), ylim=c(0.2, 0.7))
}

testPlot2 <- function() {
  td <- zeros(2,6)
  td[1,] <- c(.3, .6, .8, .1, .6, 1)
  td[2,] <- c(.3, .6, .8, .25, .6, .85)
  N <- 30  
  nd <- list(zeros(N,6),zeros(N,6))
  for (i in 1:2) {
    for (j in 1:N) {
      nd[[i]][j,] <- td[i,] + 0.5*(runif(6,-1,1)) 
    }
  }

  staPLOT(nd)
}

