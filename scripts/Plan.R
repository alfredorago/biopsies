### Plan

plan = drake_plan(
  
  ## Import and stamilize gene expression from experiment
  # Create varance stabilized gene count table
  gcnt.vst = gcnt[,grepl(pattern = "_J_", x = colnames(gcnt))] %>%
    round(.) %>%
    DESeq2::varianceStabilizingTransformation(.),
  gvst.df = reshape::melt(gcnt.vst) %>% 
    mutate(.,
           lib = str_replace(X2, "_(L|P).*", ""),
           cond = str_replace_all(str_extract(X2, "_(L|P)_"), "_", "")
    ),
  
  ## Create lookup table for gene symbols
  # Get lookup table for gene synonyms: necessary because Smillie data includes gene aliases
  aliasSymbol = dbGetQuery(conn = org.Hs.eg_dbconn(), 
                           statement = 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
  ),
  
  
  # Get list of gene symbol names and aliases for all genes in the single cell dataset
  # Necessary because several of the symbols are aliases and we need their real names too
  geneSymbols.full = 
    unique(c(aliasSymbol[which(aliasSymbol$symbol%in%rownames(smiLogCpmCorrected)),"symbol"],
                              rownames(smiLogCpmCorrected))) %>%
    checkGeneSymbols(.) %>%
    .$Suggested.Symbol,
  
  # query for ENSEMBL IDs
  ensemblIDs = getBM(attributes = list("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version", 'hgnc_symbol'), 
                     filters = "hgnc_symbol", 
                     values = geneSymbols.full,
                     mart = useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')),
  
  ## dtangle workflow
  # Import reference single-cell expression profiles
  # Rename gene IDs in expression table, and convert to matrix
  singleCellEnsembl =   
    checkGeneSymbols(rownames(smiLogCpmCorrected))[,"Suggested.Symbol"] %>% # losing 187 genes
    match(x = ., table = ensemblIDs$hgnc_symbol) %>%
    ensemblIDs[.,'hgnc_symbol'] %>%
    mutate(smiLogCpmCorrected, ensemblID = .)  %>%
    na.exclude(.) %>% # Remove genes that have no corresponding ENSEMBL ID
    .[-which(duplicated(.$"ensemblID")),] %>% # Remove genes that have more than one ENSEMBL ID
    remove_rownames(.) %>%
    column_to_rownames(., var = "ensemblID"),
  
  # Merge experiment and reference datasets using ensembl IDs 
  exprAndReference = merge(x = gcnt.vst, y = singleCellEnsembl, 
                           by = 0, sort = F,
                           all.x = F, all.y = F, 
                           suffixes = c("_test", "_refs")) %>%
    as.matrix(.) %>% 
    col2rowname(.),
  
  # Import markers from Cristoph's reanalysis of Smillie dataset, convert to Ensembl IDs and map in dataset
  markerPos = map(smillieFinalUntangledGeneList, function(x){
    names(x) %>%
      match(., ensemblIDs$external_gene_name) %>%
      ensemblIDs[. ,"ensembl_gene_id_version"] %>%
      match(., rownames(exprAndReference))
  }),
  
  # Create list of pure samples
  # Create list of names per cell types, then convert to indices using match
  refCellTypes = unique(smiCellTypes$cell_type),
  pureSamples = lapply(refCellTypes, 
                       FUN = function(x){smiCellTypes$id[which(smiCellTypes$cell_type==x)]}) %>%
    set_names(., refCellTypes) %>%
    lapply(., FUN = function(x){match(x, colnames(exprAndReference))}),
  
  # Run dTangle using singleCellEnsembl as reference
  dtangle_results = dtangle(
    Y = t(exprAndReference), 
    pure_samples = pureSamples[names(markerPos)],
    markers = markerPos,
    data_type = 'rna-seq')
  
  
  
)

