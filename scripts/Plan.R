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
    )
)
