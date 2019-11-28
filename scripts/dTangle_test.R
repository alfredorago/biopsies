### Run dtangle using full smillie markers

# Convert table to list of genes
markers <- readd(markers) %>%
  group_by(ident)

filter(.data = markers, coefD >= 2.0) %>% 
  count(.) %>% 
  full_join(count(markers), ., by = "ident") %>%
  ggplot(., aes(x = n.x, y = n.y)) + 
  geom_point() +
  geom_hline(aes(yintercept = 10))

markers.fine <- group_map(.tbl = markers, .f = ~ pull(.x, "gene")) %>% 
  set_names(., nm = group_keys(.tbl = markers)[, 1, drop=T]
  )

# Import data from smillie
load(file = "./Downloads/smillie2/smillie_downsampled_expression.Rdata")
load(file = "./Downloads/smillie2/smillie_downsampled_metadata.Rdata")

# use sql to get alias table and gene_info table (contains the symbols)
library(org.Hs.eg.db)
library(DBI)

# first open the database connection
dbCon <- org.Hs.eg_dbconn()
# write your SQL query
sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
# execute the query on the database
aliasSymbol <- dbGetQuery(dbCon, sqlQuery)
# subset to get your results
geneSymbols <- aliasSymbol[which(aliasSymbol$symbol%in%rownames(smiUmiSel)),]
# join together with genes that have no aliases
geneSymbols.full <- unique(c(geneSymbols$symbol, rownames(smiUmiSel)))

# create lookup table for symbols to ensembl GIDs
ensemblIDs <- getBM(attributes = list("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version"), 
                    filters = "external_gene_name", 
                    values = geneSymbols.full,
                    mart = useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl'))
rownames(smiUmiSel) <- 
  ensemblIDs[match(row.names(smiUmiSel), ensemblIDs$external_gene_name), "ensembl_gene_id"]



getBM(attributes = list("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version"), 
      filters = "external_gene_name", 
      values = failures,
      mart = useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl'))


## Run dtangle
cell.proportions.coarse = dtangle(
  Y = t(gcnt.vst),
  markers = markerPos.coarse,
  pure_samples = NULL,
  references = NULL,
  data_type = 'rna-seq'
)


