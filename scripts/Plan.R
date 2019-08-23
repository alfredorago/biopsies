### Plan

plan = drake_plan(
  
  ## File preprocessing
  
  gcnt.vst = gcnt[,grepl(pattern = "_J_", x = colnames(gcnt))] %>%
    round(.) %>%
    DESeq2::varianceStabilizingTransformation(.),
  gvst.df = reshape::melt(gcnt.vst) %>% 
    mutate(.,
           lib = str_replace(X2, "_(L|P).*", ""),
           cond = str_replace_all(str_extract(X2, "_(L|P)_"), "_", "")
    ),
  finalGcnt = gcnt.vst[which(row.names(readd(gcnt.vst))%in%row.names(readd(finalfc))),],
  fcvst.df = reshape::cast(gvst.df, X1+lib~cond, value="value") %>% 
    mutate(.,
           Lcnt = 2^L,
           Pcnt = 2^P,
           log2abd = log2(Lcnt + Pcnt),
           log2fc = L - P,
           part = str_replace(lib, ".*_", "")
    ),
  # removing intermediate data.frame finaldiff, which is the filtered version of fcvst.df only
  # creates tables with fold change differences
  finalfc = filter(fcvst.df, abs(log2fc) >= 1.4) %>%
    reshape::cast(., X1~lib, value="log2fc") %>%
    rows2cols(.,"X1"),
  finalfc.t  = t(data.frame(finalfc)),
  totalfc = reshape::cast(fcvst.df, X1~lib, value="log2fc") %>%
    rows2cols("X1"),
  
  # dtangle workflow
# Import and then subset marker gene set from Smillie
 xlPath = here("Downloads/smillie/SM/1-s2.0-S0092867419307329-mmc2.xlsx"),
# file_in(xlPath),

# file_out(here("./output/SmillieMarkers.csv")),
 markers = excel_sheets(path = xlPath) %>%
   set_names(.) %>%
   .[1:3] %>%
   map_dfr(.x = ., .f = ~ read_excel(path = xlPath, sheet = .x, col_types = c(
     ident = "text",
     gene = "text",
     coefD	= "numeric",
     pvalD	= "numeric",	
     padjD	= "numeric",
     coefC	= "numeric",
     pvalC	= "numeric",
     padjC	= "numeric",
     mastfc = "numeric",
     alpha	= "numeric",
     ref_alpha	= "numeric",
     mu	= "numeric",
     ref_mu	= "numeric",
     mean	= "numeric",
     ref_mean	= "numeric",
     log2fc	= "numeric",
     spec_h = "logical",
     spec_d = "logical"
   )))  %>% 
  group_by(ident) %>%
  write_csv(., path = here("./output/SmillieMarkers.csv")),

markers.filtered = filter(.data = markers, spec_h == TRUE),

markerIDs.fine = group_map(.tbl = markers, .f = ~ pull(.x, "gene")) %>% 
  set_names(., nm = group_keys(.tbl = markers)[, 1, drop=T]
  ),

# Import reference single-cell expression profiles

## Create lookup table for gene symbols
# Get lookup table for gene synonyms,
aliasSymbol = dbGetQuery(conn = org.Hs.eg_dbconn(), 
                          statement = 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
                          ),
# Get list of gene symbol names and aliases for all genes in the single cell dataset
geneSymbols.full = unique(c(aliasSymbol[which(aliasSymbol$symbol%in%rownames(smiUmiSel)),"symbol"],
                             rownames(smiUmiSel))),
# query for ENSEMBL IDs
ensemblIDs = getBM(attributes = list("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version"), 
                    filters = "external_gene_name", 
                    values = geneSymbols.full,
                    mart = useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')),
# Rename gene IDs in expression table, and convert to matrix
singleCellIndexes = list(
  ensemblIDs[match(rownames(smiUmiSel), ensemblIDs$external_gene_name),"ensembl_gene_id"],
  colnames(smiUmiSel) 
),
singleCellEnsembl = matrix(smiUmiSel, nrow = nrow(smiUmiSel), dimnames = singleCellIndexes)

# Run dTangle using singleCellEnsembl as reference


)
