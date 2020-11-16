#/*==========================================================================#*/
#' ## Extended Data Figure 8C
#/*==========================================================================#*/
# SOM promoter class expression correlation

## Housekeeping TPM correlation
# load all necesssary data
  # load soms consensus clusters
    consensus.clusters.som_promOnly.anno.grl <- readRDS("intermediate_data/consensus_clusters_som_promOnly_anno_Damir.9classes_grl.RDS")

# select housekeeping/ubiquitously expressed class - som1
  housekeep_som9.gr <- consensus.clusters.som_promOnly.anno.grl[[1]]

# load data
  domTSS_PGCs_oocyte_embryo_anno.grl <- readRDS(file ="../all_reps_PN6_2cell/intermediate_data/domTSS_PGCs_oocyte_embryo_anno_grl.RDS")

# exclude weird chromosomes
  domTSS_PGCs_oocyte_embryo_anno.grl <- lapply(domTSS_PGCs_oocyte_embryo_anno.grl, function(x){
    gr <- dropSeqlevels(x, "chrM", pruning.mode = "coarse")
  })

# select dominant TSSs in eatch sample overlappying the housekeeping SOM class of tag clusters
domTSS_housekeep_som9.grl <- lapply(domTSS_PGCs_oocyte_embryo_anno.grl, function(x) {
  overlap <- findOverlapPairs(x, housekeep_som9.gr)
  overlap.gr <- overlap@first
  
  # attach annotation
  mcols(overlap.gr) <-cbind(mcols(overlap.gr), "annotation" = overlap@second$annotation) 
  return(overlap.gr)
})

# it may be beneficial to select only 1 domTSS per consensus cluster - highest expressed but I am not sure now... so I continue with all
  sapply(domTSS_housekeep_som9.grl, length) # similar numbers of domTSS in each sample

# convert to TCs
  TC_housekeep_som9.grl <- lapply(domTSS_housekeep_som9.grl, function(x) {
    TC.gr <- GRanges(seqnames = seqnames(x),
                     ranges = IRanges(start = x$q_0.1, 
                                      end = x$q_0.9),
                     strand = strand(x))
    mcols(TC.gr) <- mcols(x)
    seqinfo(TC.gr) <- seqinfo(x)
    return(TC.gr)
  })

#----- E11.5 tag clusters that overlap with E16.5F
  E11_5_E16.5F_TC.pairs <- findOverlaps(TC_housekeep_som9.grl[["E11_5"]], TC_housekeep_som9.grl[["E16_5F"]])
  E11_5_E16.5F_TC.pairs.gr <- TC_housekeep_som9.grl[["E11_5"]][queryHits(E11_5_E16.5F_TC.pairs)] 
  E11_5_E16.5F_TC.pairs.gr$subjectHits <- subjectHits(E11_5_E16.5F_TC.pairs)
  E16_5F_E11.5_TC.pairs.gr <- TC_housekeep_som9.grl[["E16_5F"]][subjectHits(E11_5_E16.5F_TC.pairs)]

# correlation plots - tpm correlations - 11.5 with 16.5F and P6 pairwise
# prepare dataframe
  E11_5_E16_5F_TPM.df <- data.frame("E11_5_TC_tpm" = E11_5_E16.5F_TC.pairs.gr$tpm, 
                                    "E16_5_TC_tpm" = E16_5F_E11.5_TC.pairs.gr$tpm)

# plot scatterplot
  Rcorr <- round(cor(E11_5_E16_5F_TPM.df$E11_5_TC_tpm, E11_5_E16_5F_TPM.df$E16_5_TC_tpm, "method" = "spearman"), digits = 2)

  pA <- ggplot(E11_5_E16_5F_TPM.df) +
    geom_point(shape = 21, aes(x =  log2(E11_5_TC_tpm + 1), y = log2(E16_5_TC_tpm + 1)), size = 2, alpha = 0.1, shape = 21, fill = "black") +
    theme_light() +
    theme(text = element_text(size = 20, colour = "black"), 
          axis.text.x = element_text(size = 18, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black"),
          legend.title = element_blank()) +
    guides(colour = guide_legend(reverse=T)) +
    scale_y_continuous(name = "log2(E16.5F tpm + 1)", limits = c(0, 15), seq(0, 15, 5)) +
    scale_x_continuous(name = "log2(E11.5 tpm +1)", limits = c(0, 15), seq(0, 15, by = 5)) +
    geom_abline(x = 0:15, y = 0:15, lty = "dashed", colour = "red", size = 0.8) +
    annotate(geom = "text", x = 1, y = 15, label = paste0("R = ", Rcorr), color = "black", size = 6)
  
  pdf("Figures/TPM_corr_11_5_E16_5F.pdf", width = 6, height = 6)
    print(pA)
  dev.off() 


#-----E11.5 tag clusters that overlap with P6
  E11_5_P6_TC.pairs <- findOverlaps(TC_housekeep_som9.grl[["E11_5"]], TC_housekeep_som9.grl[["PN6"]])
  E11_5_P6_TC.pairs.gr <- TC_housekeep_som9.grl[["E11_5"]][queryHits(E11_5_P6_TC.pairs)] 
  E11_5_P6_TC.pairs.gr$subjectHits <- subjectHits(E11_5_P6_TC.pairs)
  P6_E11.5_TC.pairs.gr <- TC_housekeep_som9.grl[["PN6"]][subjectHits(E11_5_P6_TC.pairs)]

# correlation plots - tpm correlations - 11.5 with 16.5F and P6 pairwise
  # prepare dataframe
    E11_5_P6_TPM.df <- data.frame("E11_5_TC_tpm" = E11_5_P6_TC.pairs.gr$tpm, 
                                  "P6_TC_tpm" = P6_E11.5_TC.pairs.gr$tpm)

# plot scatterplot
  Rcorr <- round(cor(E11_5_P6_TPM.df$E11_5_TC_tpm, E11_5_P6_TPM.df$P6_TC_tpm, "method" = "spearman"), digits = 2)

  pB <- ggplot(E11_5_P6_TPM.df) +
    geom_point(shape = 21, aes(x =  log2(E11_5_TC_tpm + 1), y = log2(P6_TC_tpm + 1)), size = 2, alpha = 0.1, shape = 21, fill = "black") +
    theme_light() +
    theme(text = element_text(size = 20, colour = "black"), 
          axis.text.x = element_text(size = 18, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black"),
          legend.title = element_blank()) +
    guides(colour = guide_legend(reverse=T)) +
    scale_y_continuous(name = "log2(P6 tpm + 1)", limits = c(0, 15), seq(0, 15, 5)) +
    scale_x_continuous(name = "log2(E11.5 tpm +1)", limits = c(0, 15), seq(0, 15, by = 5)) +
    geom_abline(x = 0:15, y = 0:15, lty = "dashed", colour = "red", size = 0.8) +
    annotate(geom = "text", x = 1, y = 15, label = paste0("R = ", Rcorr), color = "black", size = 6)
  
  pdf("Figures/TPM_corr_11_5_P6.pdf", width = 6, height = 6)
    print(pB)
  dev.off() 


#------E16.5F tag clusters that overlap with P6
  E16_5F_P6_TC.pairs <- findOverlaps(TC_housekeep_som9.grl[["E16_5F"]], TC_housekeep_som9.grl[["PN6"]])
  E16_5F_P6_TC.pairs.gr <- TC_housekeep_som9.grl[["E16_5F"]][queryHits(E16_5F_P6_TC.pairs)] 
  E16_5F_P6_TC.pairs.gr$subjectHits <- subjectHits(E16_5F_P6_TC.pairs)
  P6_E16_5F_TC.pairs.gr <- TC_housekeep_som9.grl[["PN6"]][subjectHits(E16_5F_P6_TC.pairs)]

# correlation plots - tpm correlations - 11.5 with 16.5F and P6 pairwise
# prepare dataframe
  E16_5F_P6_TPM.df <- data.frame("E16_5F_TC_tpm" = E16_5F_P6_TC.pairs.gr$tpm, 
                                 "P6_TC_tpm" = P6_E16_5F_TC.pairs.gr$tpm)

# plot scatterplot
  Rcorr <- round(cor(E16_5F_P6_TPM.df$E16_5F_TC_tpm, E16_5F_P6_TPM.df$P6_TC_tpm, "method" = "spearman"), digits = 2)

  pC<- ggplot(E16_5F_P6_TPM.df) +
    geom_point(shape = 21, aes(x =  log2(E16_5F_TC_tpm + 1), y = log2(P6_TC_tpm + 1)), size = 2, alpha = 0.1, shape = 21, fill = "black") +
    theme_light() +
    theme(text = element_text(size = 20, colour = "black"), 
          axis.text.x = element_text(size = 18, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black"),
          legend.title = element_blank()) +
    guides(colour = guide_legend(reverse=T)) +
    scale_y_continuous(name = "log2(P6 tpm + 1)", limits = c(0, 15), seq(0, 15, 5)) +
    scale_x_continuous(name = "log2(E16.5F tpm +1)", limits = c(0, 15), seq(0, 15, by = 5)) +
    geom_abline(x = 0:15, y = 0:15, lty = "dashed", colour = "red", size = 0.8) +
    annotate(geom = "text", x = 1, y = 15, label = paste0("R = ", Rcorr), color = "black", size = 6)

  pdf("Figures/TPM_corr_16_5F_P6.pdf", width = 6, height = 6)
    print(pC)
  dev.off() 

  library(cowplot)
  pdf("Figures/TPM_corr_all.pdf", width = 18, height = 6)
    plot_grid(pA, pB, pC, ncol = 3)
  dev.off()


#------mESC tag clusters that overlap with P14
  mESC_P14_TC.pairs <- findOverlaps(TC_housekeep_som9.grl[["E14_mESC"]], TC_housekeep_som9.grl[["PN14"]])
  mESC_P14_TC.pairs.gr <- TC_housekeep_som9.grl[["E14_mESC"]][queryHits(mESC_P14_TC.pairs)] 
  mESC_P14_TC.pairs.gr$subjectHits <- subjectHits(mESC_P14_TC.pairs)
  P14_mESC.pairs.gr <- TC_housekeep_som9.grl[["PN14"]][subjectHits(mESC_P14_TC.pairs)]

# correlation plots - tpm correlations -mESC and PN14 and MII pairwise
  # prepare dataframe
    mESC_P14_TPM.df <- data.frame("mESC_TC_tpm" = mESC_P14_TC.pairs.gr$tpm, 
                                  "P14_TC_tpm" = P14_mESC.pairs.gr$tpm)

# plot scatterplot
  Rcorr <- round(cor(mESC_P14_TPM.df$mESC_TC_tpm, mESC_P14_TPM.df$P14_TC_tpm, "method" = "spearman"), digits = 2)

  pD<- ggplot(mESC_P14_TPM.df) +
    geom_point(shape = 21, aes(x =  log2(mESC_TC_tpm + 1), y = log2(P14_TC_tpm + 1)), size = 2, alpha = 0.1, shape = 21, fill = "black") +
    theme_light() +
    theme(text = element_text(size = 20, colour = "black"), 
          axis.text.x = element_text(size = 18, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black"),
          legend.title = element_blank()) +
    guides(colour = guide_legend(reverse=T)) +
    scale_y_continuous(name = "log2(P14 tpm + 1)", limits = c(0, 15), seq(0, 15, 5)) +
    scale_x_continuous(name = "log2(mESC tpm +1)", limits = c(0, 15), seq(0, 15, by = 5)) +
    geom_abline(x = 0:15, y = 0:15, lty = "dashed", colour = "red", size = 0.8) +
    annotate(geom = "text", x = 1, y = 15, label = paste0("R = ", Rcorr), color = "black", size = 6)
  
  pdf("Figures/TPM_corr_mESC_P14.pdf", width = 6, height = 6)
    print(pD)
  dev.off() 


#------mESC tag clusters that overlap with MII
  mESC_MII_TC.pairs <- findOverlaps(TC_housekeep_som9.grl[["E14_mESC"]], TC_housekeep_som9.grl[["oocyte"]])
  mESC_MII_TC.pairs.gr <- TC_housekeep_som9.grl[["E14_mESC"]][queryHits(mESC_MII_TC.pairs)] 
  mESC_MII_TC.pairs.gr$subjectHits <- subjectHits(mESC_MII_TC.pairs)
  MII_mESC.pairs.gr <- TC_housekeep_som9.grl[["oocyte"]][subjectHits(mESC_MII_TC.pairs)]
  
  # correlation plots - tpm correlations -mESC and PN14 and MII pairwise
  # prepare dataframe
    mESC_MII_TPM.df <- data.frame("mESC_TC_tpm" = mESC_MII_TC.pairs.gr$tpm, 
                                  "MII_TC_tpm" = MII_mESC.pairs.gr$tpm)

# plot scatterplot
  Rcorr <- round(cor(mESC_MII_TPM.df$mESC_TC_tpm, mESC_MII_TPM.df$MII_TC_tpm, "method" = "spearman"), digits = 2)

  pE<- ggplot(mESC_MII_TPM.df) +
    geom_point(shape = 21, aes(x =  log2(mESC_TC_tpm + 1), y = log2(MII_TC_tpm + 1)), size = 2, alpha = 0.1, shape = 21, fill = "black") +
    theme_light() +
    theme(text = element_text(size = 20, colour = "black"), 
          axis.text.x = element_text(size = 18, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black"),
          legend.title = element_blank()) +
    guides(colour = guide_legend(reverse=T)) +
    scale_y_continuous(name = "log2(MII tpm + 1)", limits = c(0, 15), seq(0, 15, 5)) +
    scale_x_continuous(name = "log2(mESC tpm +1)", limits = c(0, 15), seq(0, 15, by = 5)) +
    geom_abline(x = 0:15, y = 0:15, lty = "dashed", colour = "red", size = 0.8) +
    annotate(geom = "text", x = 1, y = 15, label = paste0("R = ", Rcorr), color = "black", size = 6)

  pdf("Figures/TPM_corr_mESC_MII.pdf", width = 6, height = 6)
    print(pE)
  dev.off() 
