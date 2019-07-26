### Plan

plan = drake_plan(
  gcnt.vst = DESeq2::varianceStabilizingTransformation(round(gcnt)),
  gvst.df = reshape::melt(gcnt.vst),
  gvst.df.mut = mutate(gvst.df,
                       lib = str_replace(gvst.df$X2, "_(L|P).*", ""),
                       cond = str_replace_all(str_extract(gvst.df$X2, "_(L|P)_"), "_", "")
  ),
  fcvst.df = reshape::cast(gvst.df.mut, X1+lib~cond, value="value"),
  fcvst.df2 = mutate(fcvst.df,
                     Lcnt = 2^L,
                     Pcnt = 2^P,
                     log2abd = log2(Lcnt + Pcnt),
                     log2fc = L- P,
                     part = str_replace(lib, ".*_", "")
  )
)
