# Load base functions
library(drake)
library(here)

# Load datasets
load(here("./data/gene_grch38.rda"))
load(file = here("Downloads/smillie2/smillie_downsampled_expression.Rdata"))
load(file = here("Downloads/smillie2/smillie_downsampled_metadata.Rdata"))

# Load libraries
source(file = here("scripts", "Libraries.R"))

# Create plan
source(file = here("scripts", "Plan.R"))

# Configure
drake_config(
  plan
)
