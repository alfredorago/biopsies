### Plan

plan = drake_plan(
  
  ## Import and stamilize gene expression from experiment
  # Filter only J samples, remove un-analyzed individuals and remove version from gene IDs
  gcnt.J = gcnt[,grepl(pattern = "_J_", x = colnames(gcnt))] %>%
    round(.) %>%
    set_rownames(str_extract(rownames(.), pattern = "^[^.]*")) %>%
    .[, !grepl(pattern =  c("R4|R18"), x = colnames(.))],
  
  gcnt.vst = gcnt.J %>%
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
  
  # Create lookup table between HGNC symbols and ensembl IDs (without version numbers)
  ensemblIDs = getBM(attributes = list("ensembl_gene_id", 'hgnc_symbol'), 
                     filters = "hgnc_symbol", 
                     values = geneSymbols.full,
                     mart = useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')),
  
  ## dtangle workflow
  # Import reference single-cell expression profiles
  # Convert gene symbols to ensembl IDs (must first update to current HGNC names)
  singleCellEnsembl =   
    checkGeneSymbols(rownames(smiLogCpmCorrected))[,"Suggested.Symbol"] %>% # losing 187 genes
    match(x = ., table = ensemblIDs$hgnc_symbol) %>%
    ensemblIDs[.,'ensembl_gene_id'] %>%
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
  
  # Check which smillie genes are missing from our dataset
  missingGenes = 
    rownames(singleCellEnsembl)%in%rownames(gcnt.vst) %>%
    not(.) %>%
    which(.) %>%
    rownames(singleCellEnsembl)[.],
  
  # Import markers from Cristoph's reanalysis of Smillie dataset, convert to Ensembl IDs and map in dataset
  # Save marker table
  markerTable = map(smillieFinalUntangledGeneList, names) %>% 
    stri_list2matrix(byrow = T) %>% 
    set_rownames(names(smillieFinalUntangledGeneList)) %>%
    as.data.frame(.) %T>%
    write.csv(., file = file_out(here("output/CellMarkers.csv")), row.names = T, na = ""),
  
  markerPosNAs = map(smillieFinalUntangledGeneList, function(x){
    names(x)  %>%
      checkGeneSymbols(.)  %>% 
      pull(.data = ., var = "Suggested.Symbol") %>%
      match(., ensemblIDs$hgnc_symbol) %>%
      ensemblIDs[. ,"ensembl_gene_id"] %>%
      match(., rownames(exprAndReference))
  }),
  
  markerPos = map(markerPosNAs, ~ na.exclude(.x)),
  
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
    pure_samples = pureSamples,
    markers = markerPos,
    summary_fn = median,
    data_type = 'rna-seq'),
  
  cellPropRaw = 
    dtangle_results$estimates %>%
    as.data.frame(.) %>%
    rownames_to_column() %>%
    gather(-1 ,key = 'ID', value = 'Proportion') %>%
    filter(., grepl("_J_", rowname)),
  
  cellProportions = 
    dtangle_results$estimates[which(grepl('_J_', rownames(dtangle_results$estimates))), ] %>%
    scale(., scale = T, center = T) %>%
    as.data.frame(.) %>% 
    rownames_to_column(., var = "sampleID") %>%
    gather(., key = 'cellType', value = 'proportion', -sampleID),
  
  tidyProps = 
    mutate(cellProportions,
           subj = str_extract(sampleID, pattern = "[0-9]{1,2}") %>%
             as.factor(),
           cond = str_replace_all(str_extract(sampleID, "_(L|P)_"), "_", "") %>%
             factor(levels = c('P', 'L')),
           week = str_extract(sampleID, pattern = "W[:digit:]") %>%
             factor()
    ),
  
  metadata = 
    read_tsv(file = here("Downloads/Clinical_characteristics_Herlevstudy.txt")) %>%
    dplyr::select(., contains("Screening"), contains("Treatment"), contains("Group")) %>%
    gather(., key = "Visit", value = "Treatment", ... = contains("Treatment")) %>%
    mutate(., 
           Visit = str_extract(string = Visit, pattern = "[0-9]") %>% 
             paste0("W", .) %>% 
             as.factor(.),
           Treatment = str_extract(string = Treatment, pattern = "^[A,P]") %>%
             str_replace(., "A", "L") %>%
             as.factor(.)) %>%
    set_colnames(., c('subj', 'respond', 'week', 'cond')),
  
  metaProps = 
    merge(metadata, tidyProps)

)

