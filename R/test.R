source('STA.R')

cmrData <- readMat('../data/nakabayashi.mat')
cmrData <- cmrData$data
a <- CMRfits(1, cmrData)
print(a$p)

