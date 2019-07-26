### Plan

plan = drake_plan(
  gcnt.vst = DESeq2::varianceStabilizingTransformation(round(gcnt)),
  gvst.df = reshape::melt(gcnt.vst) %>% 
    mutate(.,
           lib = str_replace(X2, "_(L|P).*", ""),
           cond = str_replace_all(str_extract(X2, "_(L|P)_"), "_", "")
    ),
  fcvst.df = reshape::cast(gvst.df, X1+lib~cond, value="value") %>% 
    mutate(.,
           Lcnt = 2^L,
           Pcnt = 2^P,
           log2abd = log2(Lcnt + Pcnt),
           log2fc = L - P,
           part = str_replace(lib, ".*_", "")
    ),
  # removing intermediate data.frame finaldiff, which is the filtered version of fcvst.df only
  finalfc = filter(fcvst.df, abs(log2fc) >= 1.4) %>%
    reshape::cast(., X1~lib, value="log2fc") %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames("X1"),
  finalfc.t  = t(data.frame(finalfc))
)
