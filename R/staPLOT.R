source('sta.R')

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
