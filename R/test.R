source('STA.R')

cmrData <- readMat('../data/nakabayashi.mat')
cmrData <- cmrData$data
CMRfits(2, cmrData)
