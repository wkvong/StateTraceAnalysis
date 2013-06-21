source('STA.R')

x <- matrix(scan('../data/x.dat'), ncol = 5, byrow = TRUE)

output <- staSTATS(x)
print(output)

## cmrData <- readMat('../data/nakabayashi.mat')
## cmrData <- cmrData$data
## a <- CMRfits(1, cmrData)
## print(a$p)

