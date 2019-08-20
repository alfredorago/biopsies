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
  # specify which genes are markers for tissues of interest (after Smillie 2019)
  markerNames.coarse = list(epithelial = c("EPCAM", "KRT8", "KRT18"),
                            stromal = c("COL1A1", "COL1A2", "COL6A1", "COL6A2", "VWF", "PLVAP", "CDH5", "S100B"),
                            immune = c("CD52", "CD2", "CD3D", "CD3G", "CD3E", "CD79A", "CD79B", "CD14", "CD16", "CD68", "CD83", "CSF1R", "FCER1G")
  ),
  markerTable.coarse = lapply(markerNames.coarse, function(x){
    getBM(attributes = list("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version"), 
          filters = "external_gene_name", 
          values = x,
          mart = useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl'))
  }),
 markerIDs.coarse = map(.x = markerTable.coarse, .f = function(x){x[,"ensembl_gene_id_version"]}), 
 markerPos.coarse = map(.x = markerIDs.coarse, .f = function(x){
   which(rownames(gcnt.vst)%in%x)
 }),
 # Using mock set of reference samples sampled at random from main dataset
 references.coarse = gcnt.vst[,sample(1:ncol(gcnt.vst), length(markerPos.coarse))] %>% 
   t(.),
 
 # Run dtangle
 dtOut.coarse = dtangle(Y=t(gcnt.vst), 
                        references = references.coarse, 
                        markers = markerPos.coarse, 
                        data_type = 'rna-seq'),
 
 # Import cell metadata from smillie dataset
 smillie.meta = read_tsv(file = "./Downloads/smillie/all.meta2.txt", 
                     col_types = 'cfnnffff', skip = 2, 
                     col_names = c("NAME", "Cluster", "nGene", "nUMI", "Subject", "Health", "Location", "Sample"))
 
 
)