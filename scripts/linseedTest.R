library(drake)
library(linseed)
library(ggplot2)
library(here)
library(magrittr)
source(file = "./scripts/unbundled.R")

pdf()

lo =readd(totalfc) %>%
 LinseedObject$new(.)
lo$calculatePairwiseLinearity()
lo$calculateSpearmanCorrelation()
lo$calculateSignificanceLevel(100)
lo$significancePlot(0.01)
lo$svdPlot()
lo$setCellTypeNumber(6)
lo$project('full') # projecting full dataset

lo$projectionPlot()



lo$smartSearchCorners(error="norm")
# lets select 100 genes closest to the simplex corners
lo$selectGenes(100)
lo$tsnePlot()


dev.off()