## Drake plan for cell deconvolution analysis

### Plan

plan = drake_plan(
  cellMarkers = data.table::fread(here("Downloads", "Human_cell_markers_Zhang2018.txt"), drop = c(1,11:15), sep2 = "auto") %>%
    .[tissueType %in% c("Blood", "Blood Vessel", "Colorectum", "Stomach", "Lymph node", "Colon", "Large Intestine", "Small Intestine", "Epithelium", "Gastrointestinal tract", "Intestine", "Jejunum", "Small intestinal crypt") & cancerType == "Normal"],
  cellReference = fread(file = here("Downloads", "hg38.cage_peak_phase1and2combined_fair_counts_ann.osc.txt.gz.extract.tsv"))
)