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
    rows2cols("X1")
)
