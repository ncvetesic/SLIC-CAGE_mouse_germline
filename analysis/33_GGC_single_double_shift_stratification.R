#/*==========================================================================#*/
#' ## Extended Data Figure 14A
#/*==========================================================================#*/
# stratification of GGC (E16.5F) vs mESC shifting promoters into single and double shift

# Overlap of E16.5F vs mESC shifting and E16.5F vs oocyte shifting - do late GGC shifts shift to maternal?
# import shifting ranges
  PGC_oocyte_embryo_vs_mESC_shift.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/PGC_oocyte_embryo_vs_mESC_shift_grl.RDS")
  PGC_oocyte_embryo_vs_E16_5.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/PGC_oocyte_embryo_vs_E16_5_grl.RDS")

# select 16_5F shifts
  E16_5F_vs_mESC_shifts.gr <- PGC_oocyte_embryo_vs_mESC_shift.grl[["E16_5F"]]
  PN6_vs_E16_5F_shifts.gr <- PGC_oocyte_embryo_vs_E16_5.grl[["PN6"]]

# overlap shifts oocyte vs E16_5F and E16_5F vs mESC - this will get which shift again - but wont know if they shift to oocyte or somatic code
# in double_shift_E16_5F_vs_mESC, X group == E16_5F and Y group == mESC
  double_shift_E16_5F_vs_mESC.gr <- subsetByOverlaps(E16_5F_vs_mESC_shifts.gr, PN6_vs_E16_5F_shifts.gr, maxgap = 0)
  double_shift_E16_5F_vs_mESC.gr[order(double_shift_E16_5F_vs_mESC.gr$shifting.score, decreasing = T)][15:25]

# in double_shift_PN6_vs_E16_5F, X group == PN6 and Y group == mESC
  double_shift_PN6_vs_E16_5F.gr <- subsetByOverlaps(PN6_vs_E16_5F_shifts.gr, E16_5F_vs_mESC_shifts.gr, maxgap = 0)
  double_shift_PN6_vs_E16_5F.gr[order(double_shift_PN6_vs_E16_5F.gr$shifting.score, decreasing = T)][15:25]

# to ensure distinguish if they shift to somatic or oocyte TSSs i should overlap with  PN6_vs_mESC shifts
# select PN6_vs_mESC shifts
  PN6_vs_mESC_shifts.gr <- PGC_oocyte_embryo_vs_mESC_shift.grl[["PN6"]]

# overlap with PN6 vs mESC shifts - shitfing GGCs that shift to maternal - max gap can be 10bp because these are consensus clusters
  double_shift_E16_5F_vs_mESC_maternal.gr <- subsetByOverlaps(double_shift_E16_5F_vs_mESC.gr, PN6_vs_mESC_shifts.gr, maxgap = 0)
  double_shift_E16_5F_vs_mESC_maternal.gr[order(double_shift_E16_5F_vs_mESC_maternal.gr$shifting.score, decreasing = T)][1:10]

# shifting GGCs that revert to somatic
  overlap.pairs <- findOverlaps(double_shift_E16_5F_vs_mESC.gr, PN6_vs_mESC_shifts.gr, maxgap = 0)
  query.hits <- queryHits(overlap.pairs)
  double_shift_E16_5F_vs_mESC_somatic.gr <- double_shift_E16_5F_vs_mESC.gr[-query.hits]
  double_shift_E16_5F_vs_mESC_somatic.gr[order(double_shift_E16_5F_vs_mESC_somatic.gr$shifting.score, decreasing = T)][10:20]
  rm(overlap.pairs)
  rm(query.hits)

# single shift - is they dont shift in P6/oocyte stage
  doubleshift.pairs <- findOverlaps(E16_5F_vs_mESC_shifts.gr, PN6_vs_E16_5F_shifts.gr, maxgap = 0)
  query.hits <- queryHits(doubleshift.pairs)
  single_shift_E16_5F_vs_mESC.gr <- E16_5F_vs_mESC_shifts.gr[-query.hits]
  single_shift_E16_5F_vs_mESC.gr[order(single_shift_E16_5F_vs_mESC.gr$shifting.score, decreasing = T)][10:20]

  length(single_shift_E16_5F_vs_mESC.gr)
  length(double_shift_E16_5F_vs_mESC_maternal.gr)
  length(double_shift_E16_5F_vs_mESC_somatic.gr)

# -- all above sum up to length of total E16_5F shifts -- #
  length(E16_5F_vs_mESC_shifts.gr)

# ----- subset GGC shifting to overlap SOM housekeeping class ---- #     
# import SOM classification - som1 are housekeeping genes
  consensus.clusters.som_promOnly.anno.grl <- readRDS("intermediate_data/consensus_clusters_som_promOnly_anno_Damir.9classes_grl.RDS")
  housekeping_som1.gr <- consensus.clusters.som_promOnly.anno.grl[[1]]

# select 16_5F shifts
  E16_5F_vs_mESC_shifts.gr <- PGC_oocyte_embryo_vs_mESC_shift.grl[["E16_5F"]]
  PN6_vs_E16_5F_shifts.gr <- PGC_oocyte_embryo_vs_E16_5.grl[["PN6"]]
  PN6_vs_mESC_shifts.gr <- PGC_oocyte_embryo_vs_mESC_shift.grl[["PN6"]]

# overlap E16_5F shifting with housekeeping class first to subselect ubiquitously expressed
  E16_5F_vs_mESC_shifts_housekeep.gr <- subsetByOverlaps(E16_5F_vs_mESC_shifts.gr, housekeping_som1.gr, maxgap = 0)
  PN6_vs_E16_5F_shifts_housekeep.gr <- subsetByOverlaps(PN6_vs_E16_5F_shifts.gr, housekeping_som1.gr, maxgap = 0)
  PN6_vs_mESC_shifts_housekeep.gr <- subsetByOverlaps(PN6_vs_mESC_shifts.gr, housekeping_som1.gr, maxgap = 0)

# overlap shifts oocyte vs E16_5F and E16_5F vs mESC - this will get which shift again - but wont know if they shift to oocyte or somatic code
# in double_shift_E16_5F_vs_mESC, X group == E16_5F and Y group == mESC
  double_shift_E16_5F_vs_mESC_housekeep.gr <- subsetByOverlaps(E16_5F_vs_mESC_shifts_housekeep.gr, PN6_vs_E16_5F_shifts_housekeep.gr, maxgap = 0)
  double_shift_E16_5F_vs_mESC_housekeep.gr[order(double_shift_E16_5F_vs_mESC_housekeep.gr$shifting.score, decreasing = T)][15:25]

# overlap with PN6 vs mESC shifts - shitfing GGCs that shift to maternal - max gap can be 0 because these are consensus clusters
  double_shift_E16_5F_vs_mESC_housekeep_maternal.gr <- subsetByOverlaps(double_shift_E16_5F_vs_mESC_housekeep.gr, PN6_vs_mESC_shifts_housekeep.gr, maxgap = 0)
  double_shift_E16_5F_vs_mESC_housekeep_maternal.gr[order(double_shift_E16_5F_vs_mESC_housekeep_maternal.gr$shifting.score, decreasing = T)][1:10]

# shifting GGCs that revert to somatic
  overlap.pairs <- findOverlaps(double_shift_E16_5F_vs_mESC_housekeep.gr, PN6_vs_mESC_shifts_housekeep.gr, maxgap = 0)
  query.hits <- queryHits(overlap.pairs)
  double_shift_E16_5F_vs_mESC_housekeep_somatic.gr <- double_shift_E16_5F_vs_mESC_housekeep.gr[-query.hits]
  double_shift_E16_5F_vs_mESC_housekeep_somatic.gr[order(double_shift_E16_5F_vs_mESC_housekeep_somatic.gr$shifting.score, decreasing = T)][10:20]
  rm(overlap.pairs)
  rm(query.hits)

# single shift - is they dont shift in P6/oocyte stage
  doubleshift.pairs <- findOverlaps(E16_5F_vs_mESC_shifts_housekeep.gr, PN6_vs_E16_5F_shifts_housekeep.gr, maxgap = 0)
  query.hits <- queryHits(doubleshift.pairs)
  single_shift_E16_5F_vs_mESC_housekeep.gr <- E16_5F_vs_mESC_shifts_housekeep.gr[-query.hits]
  single_shift_E16_5F_vs_mESC_housekeep.gr[order(single_shift_E16_5F_vs_mESC_housekeep.gr$shifting.score, decreasing = T)][10:20]

  length(single_shift_E16_5F_vs_mESC_housekeep.gr)
  length(double_shift_E16_5F_vs_mESC_housekeep_maternal.gr)
  length(double_shift_E16_5F_vs_mESC_housekeep_somatic.gr)

# -- all above sum up to length of total E16_5F shifts -- #
  length(E16_5F_vs_mESC_shifts_housekeep.gr)

# extract single shift that dont overlap with P6_vs_mESC
  singleshift.pairs <- findOverlaps(single_shift_E16_5F_vs_mESC_housekeep.gr, PN6_vs_mESC_shifts_housekeep.gr, maxgap = 0)
  query.hits <- queryHits(singleshift.pairs)
  single_shift_E16_5F_vs_mESC_housekeep_no.gr <- single_shift_E16_5F_vs_mESC_housekeep.gr[-query.hits]
  single_shift_E16_5F_vs_mESC_housekeep_no.gr[order(single_shift_E16_5F_vs_mESC_housekeep_no.gr$shifting.score, decreasing = T)][10:20]    

# Venn diagram of shifts - to show how many are double shifting and to which state
  library(ChIPpeakAnno)

# selection of overlaps
  sel <- list(E16_5F_vs_mESC_shifts_housekeep.gr, PN6_vs_E16_5F_shifts_housekeep.gr, PN6_vs_mESC_shifts_housekeep.gr)
  names <- c("E16_5F_vs_mESC", "P6_vs_E16_5F", "P6_vs_mESC")

  pdf("Figures/E16_5F_double_single_shifts_Venn.pdf", height = 4, width = 4)
    res <- makeVennDiagram(sel, NameOfPeaks = names, fill = c("#3FBC73FF", "#80b1d3", "#440154FF"))
  dev.off()  