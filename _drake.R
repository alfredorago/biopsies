# Load base functions
library(drake)
library(here)

# Load datasets
load(here("./data/gene_grch38.rda"))
load(file = here("Downloads/Cristoph/20_smillie_limma_corrected_balanced_subset.Rdata"))
load(file = here("Downloads/Cristoph/22_1_smillie_dTangle_deconvolution_final_markers.Rdata"))

# Load libraries
source(file = here("scripts", "Libraries.R"))

# Create plan
source(file = here("scripts", "Plan.R"))

# Configure
drake_config(
  plan
)
