### linSeed workflow
library(drake)
library(linseed)
library(magrittr)

# Load data and remove samples with missing data

gcnt <- readd(gcnt.vst)
gcnt <- gcnt[,-grep(x = colnames(gcnt), pattern = 'R18_J|R4_J')]

# Initialize linseed object
lo <- LinseedObject$new(gcnt)

# Perform linearity detection
lo$calculatePairwiseLinearity()
lo$calculateSpearmanCorrelation()
lo$calculateSignificanceLevel(iters=1000,
                              spearmanThreshold=0.1,
                              retVal=F)
lo$significancePlot(0.001)

# Select linear genes only
lo$filterDatasetByPval(0.001)

# Check components
lo$svdPlot()
varVector <- svd(lo$exp$filtered$norm)$d^2
cumsum(varVector/sum(varVector))

# Set cell number
lo$setCellTypeNumber(6)


#
lo$project("full") # projecting full dataset
lo$projectionPlot(color="filtered")


#
lo$project("filtered")
lo$smartSearchCorners(taus = 2^seq(0, -20, -1),
                      sisalIter=100,
                      dataset="filtered",
                      error="norm")

lo$deconvolveByEndpoints()
plotProportions(lo$proportions)

#
lo$selectGenes(100)
lo$tsnePlot()
