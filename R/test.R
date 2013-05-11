source('STA.R')

cmrData <- readMat('../data/nakabayashi.mat')
cmrData <- cmrData$data
CMRfits(1, cmrData)
