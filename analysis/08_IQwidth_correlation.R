#/*==========================================================================#*/
#' ## Figure 4E
#/*==========================================================================#*/
# comparison of IQ-widths of E16.5F vs mESC shifting promoters - width in E16.5F stage vs width in mESCs

library(BSgenome.Mmusculus.UCSC.mm10)

# load all tag clusters - needed because it contains sample specific IQ-widths
  tc_PGCs_oocyte_embryo_anno.grl <- readRDS( "../all_reps_PN6_2cell/intermediate_data/tc_PGCs_oocyte_embryo_anno_grl.RDS")

# load shifting tag clusters
  PGC_oocyte_embryo_vs_mESC_shift_anno.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/PGC_oocyte_embryo_vs_mESC_shift_anno_grl.RDS")


# mESC tag clusters
  mESC_shifting_TC.pairs <- findOverlaps(tc_PGCs_oocyte_embryo_anno.grl[["E14_mESC"]], PGC_oocyte_embryo_vs_mESC_shift_anno.grl[["E16_5F" ]])
  mESC_shifting_TC.pairs.gr <- tc_PGCs_oocyte_embryo_anno.grl[["E14_mESC"]][queryHits(mESC_shifting_TC.pairs)] 
  mESC_shifting_TC.pairs.gr$subjectHits <- subjectHits(mESC_shifting_TC.pairs)

library(dplyr)

# select only  tcs of max tpm value - multiple queryHits per subjectHit, I will focus on just single subjectHits with max tpm
  mESC_shifting_TC.df <- as.data.frame(mESC_shifting_TC.pairs.gr) %>% group_by(subjectHits) %>% dplyr::slice(which.max(tpm))
  mESC_shifting_TC.gr <- makeGRangesFromDataFrame(mESC_shifting_TC.df, keep.extra.columns = T)

# E16_5 tag clusters - select F only
# load tag clusters
  E16_5_shifting_TC.pairs <- findOverlaps(tc_PGCs_oocyte_embryo_anno.grl[["E16_5F"]], PGC_oocyte_embryo_vs_mESC_shift_anno.grl[["E16_5F" ]])
  E16_5_shifting_TC.gr <- tc_PGCs_oocyte_embryo_anno.grl[["E16_5F"]][queryHits(E16_5_shifting_TC.pairs)] 
  E16_5_shifting_TC.gr$subjectHits <- subjectHits(E16_5_shifting_TC.pairs)

library(dplyr)

# select only  tcs of max tpm value - multiple queryHits per subjectHit, I will focus on just single subjectHits with max tpm
  E16_5_shifting_TC.df <- as.data.frame(E16_5_shifting_TC.gr) %>% group_by(subjectHits) %>% dplyr::slice(which.max(tpm))
  E16_5_shifting_TC.gr <- makeGRangesFromDataFrame(E16_5_shifting_TC.df, keep.extra.columns = T)

# overlapping pairs of TCs mESC vs 16_5F
  mESC_vs_E16_5F.pairs <- findOverlaps(mESC_shifting_TC.gr, E16_5_shifting_TC.gr)
  mESC_vs_E16_5F.gr <- mESC_shifting_TC.gr[queryHits(mESC_vs_E16_5F.pairs)]
  mcols_E16_5F.df <- as.data.frame(E16_5_shifting_TC.gr[subjectHits(mESC_vs_E16_5F.pairs)])
  colnames(mcols_E16_5F.df) <- paste0(colnames(mcols_E16_5F.df), "_E16_5F")

  mcols(mESC_vs_E16_5F.gr) <- cbind(mcols(mESC_vs_E16_5F.gr), mcols_E16_5F.df)

# plot correlations of IQwidths
  IQ_width.df <- data.frame("mESC_IQwidth" = mESC_vs_E16_5F.gr$interquantile_width, "E16_5_IQwidth" = mESC_vs_E16_5F.gr$interquantile_width_E16_5F)

# plot scatterplot
  p <- ggplot(IQ_width.df) +
    geom_point(aes(x =  mESC_IQwidth, y = E16_5_IQwidth), size = 3, alpha = 0.3, shape = 21, fill = "black") +
    theme_light() +
    theme(text = element_text(size = 20, colour = "black"), 
          axis.text.x = element_text(size = 18, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black"),
          legend.title = element_blank()) +
    guides(colour = guide_legend(reverse=T)) +
    scale_y_continuous(name = "E16_5_IQwidth", limits = c(0, 250), seq(0, 250, 50)) +
    scale_x_continuous(name = "E14_5F_IQwidth", limits = c(0, 250), seq(0, 250, by = 50)) +
    geom_abline(x = 0:500, y = 0:500, lty = "dashed", colour = "red", size = 0.8) 

  pdf("Figures/IQ_width_corr_16_5_mESC.pdf", width = 6, height = 6)
    print(p)
  dev.off() 