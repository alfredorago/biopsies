## Load libraries for Drake
library(here)
library(data.table)
library(tidyverse)
library(stringi)
library(dtangle)
library(biomaRt)
library(readxl)
library(org.Hs.eg.db)
library(DBI)
library(magrittr)
library(HGNChelper)

## Utility functions with human-readable names

rows2cols = function(x, column) {
  tibble::remove_rownames(.data = x) %>%
    tibble::column_to_rownames(.data = ., var = column)
}

col2rowname = function(x){
  newmatrix = x[,-1]
  newmatrix = apply(newmatrix, c(1,2), as.numeric)
  newnames = x[,1]
  row.names(newmatrix) <- newnames
  newmatrix
}
