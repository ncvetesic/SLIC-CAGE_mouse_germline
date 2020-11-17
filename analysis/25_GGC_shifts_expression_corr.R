#/*==========================================================================#*/
#' ## Extended Data Figure 12G
#/*==========================================================================#*/
# Gonadal germ cell vs mESC shifting promoters - correlation of expression in GGC vs mESC stage

# E16.5F shifting promoters TPM correlations
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(dplyr)
  library(cowplot)
  
# load all tag clusters
  tc_PGCs_oocyte_embryo_anno.grl <- readRDS( "../all_reps_PN6_2cell/intermediate_data/tc_PGCs_oocyte_embryo_anno_grl.RDS")

# load shifting tag clusters
  PGC_oocyte_embryo_vs_mESC_shift_anno.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/PGC_oocyte_embryo_vs_mESC_shift_anno_grl.RDS")


# mESC tag clusters
  mESC_shifting_TC.pairs <- findOverlaps(tc_PGCs_oocyte_embryo_anno.grl[[1]], PGC_oocyte_embryo_vs_mESC_shift_anno.grl[["E16_5F" ]])
  mESC_shifting_TC.pairs.gr <- tc_PGCs_oocyte_embryo_anno.grl[[1]][queryHits(mESC_shifting_TC.pairs)] 
  mESC_shifting_TC.pairs.gr$subjectHits <- subjectHits(mESC_shifting_TC.pairs)

# select only  tcs of max tpm value - multiple queryHits per subjectHit, I will focus on just single subjectHits with max tpm
  mESC_shifting_TC.df <- as.data.frame(mESC_shifting_TC.pairs.gr) %>% group_by(subjectHits) %>% dplyr::slice(which.max(tpm))
  mESC_shifting_TC.gr <- makeGRangesFromDataFrame(mESC_shifting_TC.df, keep.extra.columns = T)

# E16_5 tag clusters - I'll select F only
  E16_5_shifting_TC.pairs <- findOverlaps(tc_PGCs_oocyte_embryo_anno.grl[["E16_5F"]], PGC_oocyte_embryo_vs_mESC_shift_anno.grl[["E16_5F" ]])
  E16_5_shifting_TC.gr <- tc_PGCs_oocyte_embryo_anno.grl[["E16_5F"]][queryHits(E16_5_shifting_TC.pairs)] 
  E16_5_shifting_TC.gr$subjectHits <- subjectHits(E16_5_shifting_TC.pairs)

# select only  tcs of max tpm value - multiple queryHits per subjectHit, I will focus on just single subjectHits with max tpm
  E16_5_shifting_TC.df <- as.data.frame(E16_5_shifting_TC.gr) %>% group_by(subjectHits) %>% dplyr::slice(which.max(tpm))
  E16_5_shifting_TC.gr <- makeGRangesFromDataFrame(E16_5_shifting_TC.df, keep.extra.columns = T)

# overlapping pairs of TCs mESC vs 16_5F
  mESC_vs_E16_5F.pairs <- findOverlaps(mESC_shifting_TC.gr, E16_5_shifting_TC.gr)
  mESC_vs_E16_5F.gr <- mESC_shifting_TC.gr[queryHits(mESC_vs_E16_5F.pairs)]
  mcols_E16_5F.df <- as.data.frame(E16_5_shifting_TC.gr[subjectHits(mESC_vs_E16_5F.pairs)])
  colnames(mcols_E16_5F.df) <- paste0(colnames(mcols_E16_5F.df), "_E16_5F")

  mcols(mESC_vs_E16_5F.gr) <- cbind(mcols(mESC_vs_E16_5F.gr), mcols_E16_5F.df)

# plot correlations of TPMs of a tag cluster
  data.df <- data.frame("mESC_tc_tpm" = mESC_vs_E16_5F.gr$tpm, "E16_5_tc_tpm" = mESC_vs_E16_5F.gr$tpm_E16_5F, 
                        "mESC_domTSS_tpm" = mESC_vs_E16_5F.gr$tpm.dominant_ctss, "E16_5_domTSS_tpm" = mESC_vs_E16_5F.gr$tpm.dominant_ctss_E16_5F,
                        "mESC_nrCTSS" = mESC_vs_E16_5F.gr$nr_ctss, "E16_5_nrCTSS" = mESC_vs_E16_5F.gr$nr_ctss_E16_5F)

# ---- plot scatterplot of TC TPM correlation
  Rcorr <- cor(data.df$mESC_tc_tpm, data.df$E16_5_tc_tpm, "method" = "spearman")

  pA <- ggplot(data.df) +
    geom_point(shape = 21, aes(x =  log2(mESC_tc_tpm + 1), y = log2(E16_5_tc_tpm +1)), size = 2, alpha = 0.2, shape = 21, fill = "black") +
    theme_light() +
    theme(text = element_text(size = 20, colour = "black"), 
          axis.text.x = element_text(size = 18, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black"),
          legend.title = element_text(size = 18, colour = "black")) +
    guides(colour = guide_legend(reverse=T)) +
    scale_y_continuous(name = "log2(E16_5_tpm +1)", limits = c(0, 20), seq(0, 20, by = 5)) +
    scale_x_continuous(name = "log2(mESC_tpm + 1)", limits = c(0, 20), seq(0, 20, by = 5)) +
    geom_abline(x = 0:20, y = 0:20, lty = "dashed", colour = "red", size = 0.8) +
    ggtitle(label = "Tag cluster expression")


# ---- plot scatterplot of dominant TSS TPMs correlations
  pB <- ggplot(data.df) +
    geom_point(shape = 21, aes(x =  log2(mESC_domTSS_tpm + 1), y = log2(E16_5_domTSS_tpm +1)), size = 2, alpha = 0.2, shape = 21, fill = "black") +
    theme_light() +
    theme(text = element_text(size = 20, colour = "black"), 
          axis.text.x = element_text(size = 18, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black"),
          legend.title = element_text(size = 18, colour = "black")) +
    guides(colour = guide_legend(reverse=T)) +
    scale_y_continuous(name = "log2(E16_5_dom_tpm +1)", limits = c(0, 20), seq(0, 20, by = 5)) +
    scale_x_continuous(name = "log2(mESC_dom_tpm + 1)", limits = c(0, 20), seq(0, 20, by = 5)) +
    geom_abline(x = 0:20, y = 0:20, lty = "dashed", colour = "red", size = 0.8) +
    ggtitle(label = "Dominant CTSS expression")

# ---- plot correlations of TPMs of a tag cluster
  pC <- ggplot(data.df) +
    geom_point(shape = 21, aes(x =  log2(mESC_nrCTSS + 1), y = log2(E16_5_nrCTSS + 1)), size = 2, alpha = 0.2, shape = 21, fill = "black") +
    theme_light() +
    theme(text = element_text(size = 20, colour = "black"), 
          axis.text.x = element_text(size = 18, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black"),
          legend.title = element_text(size = 18, colour = "black")) +
    guides(colour = guide_legend(reverse=T)) +
    geom_abline(x = 0:20, y = 0:20, lty = "dashed", colour = "red", size = 0.8) +
    ggtitle(label = "Number of CTSSs")
  
  pdf("Figures/tc_domTS_noCTSS_corr_16_5_mESC.pdf", width = 18, height = 6)
    plot_grid(pA, pB, pC, align = "h", ncol = 3)
  dev.off() 
```