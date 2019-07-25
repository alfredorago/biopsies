### Precluster script
# Import data and replicate PCA/hclust (sanity check)
# Cluster based on control data only, then label by responder/non-responder
# Optional: DEseq responders vs non-responders based on ctrl data only

# Set local paths
library(here)
outdir <- here("output", "ResponderVsNonResponders_crl_only")
dir.create(outdir)
graphdir <- here("graphs", "ResponderVsNonResponders_crl_only")
dir.create(graphdir)

# Import libraries
library(data.table)
library(tidyverse)
library(ggplot2)

# Import data
load(here("data", "gene_grch38.rda"))
# ganno: gene annotation table
# gcnt == txi.gene$counts

# Subset only jejunum data

gcnt <- as.data.table(gcnt)
jSamples <- grep(x = names(gcnt), pattern = "_J_", value = T)
gcntJ <- gcnt[, ..jSamples]

# Use variance stabilizing transformation on count data
gcntJ <- as.matrix(round(gcntJ))
gvst <- DESeq2::varianceStabilizingTransformation(round(gcntJ))

# Plot PCA
pcData <- prcomp(gvst)
ggplot(data = as.data.frame(pcData$rotation), aes(x=PC1, y=PC2)) +
  geom_point()

# hclust
distData <- dist(t(gvst))
hData <- hclust(distData)
plot(hData)

### Subset only pre-treatment and repeat
gvstPlacebo <- grepl(pattern = "_P_", x = colnames(gvst)) %>%
  gvst[,.]

# PCA
pcData <- prcomp(gvstPlacebo, center = T, scale. = F)
plot(pcData)
ggplot(data = as.data.frame(pcData$rotation), aes(x=PC1, y=PC2)) +
  geom_point()

# hclust
distData <- dist(t(gvstPlacebo))
hData <- hclust(distData)
plot(hData)