# Load base functions
library(drake)
library(here)

# Load datasets
load(here("./data/gene_grch38.rda"))

# Load libraries
source(file = here("scripts", "Libraries.R"))

# Create plan
source(file = here("scripts", "Plan.R"))

# Configure
drake_config(
  plan
)
