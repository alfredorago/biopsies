## Load libraries for Drake
library(here)
library(data.table)
library(tidyverse)
library(dtangle)
library(biomaRt)

## Utility functions with human-readable names

rows2cols = function(x, column) {
  tibble::remove_rownames(.data = x) %>%
    tibble::column_to_rownames(.data = ., var = column)
}
