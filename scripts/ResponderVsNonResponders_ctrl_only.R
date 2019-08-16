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
library(drake)

# Import data
function(x) load(here("data", "gene_grch38.rda"))
# ganno: gene annotation table
# gcnt == txi.gene$counts

# # Subset only jejunum data
# gcnt <- as.data.table(gcnt)
# jSamples <- grep(x = names(gcnt), pattern = "_J_", value = T)
# gcntJ <- gcnt[, ..jSamples]

#### try using DESEQ to correct fold changes (from Yun)
gcnt.vst <- DESeq2::varianceStabilizingTransformation(round(gcnt))
gvst.df <- melt(gcnt.vst)
gvst.df$lib <- str_replace(gvst.df$Var2, "_(L|P).*", "")
gvst.df$cond <- str_replace_all(str_extract(gvst.df$Var2, "_(L|P)_"), "_", "")

fcvst.df <- reshape::cast(gvst.df, Var1+lib~cond, value="value")

fcvst.df$Lcnt <- 2^fcvst.df$L
fcvst.df$Pcnt <- 2^fcvst.df$P
fcvst.df$log2abd <- log2(rowSums(fcvst.df[, c("Lcnt", "Pcnt")]))
fcvst.df$log2fc <- fcvst.df$L- fcvst.df$P

fcvst.df$part <- str_replace(fcvst.df$lib, ".*_", "")

qbase <- ggplot(data=fcvst.df)
qhex <- geom_hex(aes( x=log2abd, y=log2fc),
                 binwidth=c((range(fcfinal.df$log2abd)[2] -  range(fcfinal.df$log2abd)[1] )/250,
                            (range(fcfinal.df$log2fc)[2] -  range(fcfinal.df$log2fc)[1] )/250))

## using 2^1.4 FC as a cutoff for subpopulational studies
finaldiff <- subset(fcvst.df, abs(log2fc)>=1.4)
finalfc <- fcvst.df[chr(fcvst.df$X1) %in% chr(finaldiff$X1), ]
finalfc <- cast(finalfc, X1~lib, value="log2fc")
finalfc <- col2rowname(finalfc)
finalfc.t  <- t(data.frame(finalfc))

save(finaldiff, finalfc, finalfc.t, file="fc_subpopulation_log2fc1_4.rda")
