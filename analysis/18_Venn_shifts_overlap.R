#/*==========================================================================#*/
#' ## Extended Data Figure 4E, 12E
#/*==========================================================================#*/
# Overlap of shifting promoters sets - Venn diagrams

#### E16.5F-E16.5 P6 4-cell shifts Venn diagram of shifting overlaps and overlap granges
library(ChIPpeakAnno)

# read in data
  PGC_oocyte_embryo_vs_mESC_shift.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/PGC_oocyte_embryo_vs_mESC_shift_grl.RDS")

# selection of shifting
  sel <- PGC_oocyte_embryo_vs_mESC_shift.grl[c("E16_5F",  "E16_5M", "PN6", "S4_cell")]

  pdf("Figures/shifting_prom_oocyte_GGC_Venn.pdf", height = 4, width = 4)
    res <- makeVennDiagram(sel, NameOfPeaks = names(sel), fill = c("#3FBC73FF", "#31B57BFF" , "#440154FF", "#2E6E8EFF"))
  dev.off()

# GRanges common shifting promoters
  GGC_oocyte_embryo_shifting_overlaps.gr <- Reduce(intersect, sel)


