### Import datasets from yun's R format
library(here)
outdir <- here("output", "ImportData")
dir.create(outdir)

# Import data
load(here("data", "gene_grch38.rda"))

### Save as csv
# ganno: gene annotation table
write.csv(x = ganno, file = file.path(outdir, "GeneAnnotation.csv"), row.names = F)

# gcnt == txi.gene$counts
# Contains gene counts
write.csv(x = gcnt, file = file.path(outdir, "GeneCounts.csv"), row.names = F)
