#/*==========================================================================#*/
#' ## Extended Data Figure 4D, 12D
#/*==========================================================================#*/
# tpm distribution in shifting promoters

# import data
  library(BSgenome.Mmusculus.UCSC.mm10)
  domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/domTSS_PGC_oocyte_embryo_vs_mESC_shift_grl.RDS")
  domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/domTSS_PGC_oocyte_embryo_vs_mESC_shift_shift_X_grl.RDS")
  PGC_oocyte_embryo_vs_mESC_shift_anno.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/PGC_oocyte_embryo_vs_mESC_shift_anno_grl.RDS")

  domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl <- lapply(domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl, function(x){
    gr <- dropSeqlevels(x, "chrM", pruning.mode = "coarse")
  })

  domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl <- lapply(domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl, function(x){
    gr <- dropSeqlevels(x, "chrM", pruning.mode = "coarse")
  })

  PGC_oocyte_embryo_vs_mESC_shift_anno.grl <- lapply(PGC_oocyte_embryo_vs_mESC_shift_anno.grl, function(x){
    gr <- dropSeqlevels(x, "chrM", pruning.mode = "coarse")
  })

# load annotated individual all TCs - I need cluster tpm info
  tc_PGCs_oocyte_embryo_anno.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/tc_PGCs_oocyte_embryo_anno_grl.RDS")

  tc_PGCs_oocyte_embryo_anno.grl <- lapply(tc_PGCs_oocyte_embryo_anno.grl, function(x){
    gr <- dropSeqlevels(x, "chrM", pruning.mode = "coarse")
  })

# use just the shifting consensus cluster
  # extract cluster tpm values - groupX = x, groupY = "E14_mESC"
  shifting_tpm.dfl <- lapply(1:length(PGC_oocyte_embryo_vs_mESC_shift_anno.grl), function(x) {
      df <- data.frame("TC_tpm" = PGC_oocyte_embryo_vs_mESC_shift_anno.grl[[x]]$groupX.tpm,
                       "sample" = rep(names(PGC_oocyte_embryo_vs_mESC_shift_anno.grl)[[x]], 
                                    times = length(PGC_oocyte_embryo_vs_mESC_shift_anno.grl[[x]])))
      return(df)
  })


# calculate tpm distribution in all tag clusters in each sample
  # extract all tag cluster tpm values - exclude mESC as all comparisons are against mESC
    all_tc_tpm.dfl <- lapply(1:length(tc_PGCs_oocyte_embryo_anno.grl[-1]), function(x) {
      df <- data.frame("TC_tpm"= tc_PGCs_oocyte_embryo_anno.grl[-1][[x]]$tpm,
                       "sample" = rep(names(tc_PGCs_oocyte_embryo_anno.grl[-1])[[x]], 
                                      times = length(tc_PGCs_oocyte_embryo_anno.grl[-1][[x]])))
      df$sample <- paste0(df$sample,"_all_TCs")
      return(df)
    })

# flatten to a dataframe
  all_combined_tpm.df <- rbind(do.call(rbind, shifting_tpm.dfl), do.call(rbind, all_tc_tpm.dfl))

#-----select 16.5F/M, PN6, oocyte, 2-cell, 4-cell and all tc
  all_combined_tpm.df <- all_combined_tpm.df[all_combined_tpm.df$sample == "E16_5F" | all_combined_tpm.df$sample == "E16_5M" | all_combined_tpm.df$sample == "PN6" | all_combined_tpm.df$sample == "PN14" | all_combined_tpm.df$sample == "oocyte" | all_combined_tpm.df$sample == "S2_cell" | all_combined_tpm.df$sample == "S4_cell" | all_combined_tpm.df$sample == "E16_5F_all_TCs" | all_combined_tpm.df$sample == "E16_5M_all_TCs" | all_combined_tpm.df$sample == "PN6_all_TCs" | all_combined_tpm.df$sample == "PN14_all_TCs" | all_combined_tpm.df$sample == "oocyte_all_TCs" | all_combined_tpm.df$sample == "S2_cell_all_TCs" | all_combined_tpm.df$sample == "S4_cell_all_TCs", ]

# drop unused sequence levels
  all_combined_tpm.df$sample <- droplevels(all_combined_tpm.df$sample)

# set levels for plotting
  all_combined_tpm.df$sample <- factor(all_combined_tpm.df$sample, levels = c("E16_5F", "E16_5F_all_TCs", "E16_5M", "E16_5M_all_TCs", "PN6", "PN6_all_TCs", "PN14", "PN14_all_TCs", "oocyte", "oocyte_all_TCs", "S2_cell", "S2_cell_all_TCs", "S4_cell", "S4_cell_all_TCs"))

  col <- col_scheme[c("E16.5F", "E16.5M", "P6", "P14", "MII", "2-cell", "4-cell")]
  samples <- levels(all_combined_tpm.df$sample)

# ggplot 
  p <- ggplot(all_combined_tpm.df, aes(y = log2(TC_tpm +1), x = sample , fill = sample, alpha = 0.6)) +
    geom_boxplot(outlier.shape = NA, varwidth = T, notch = F)  +
    scale_fill_manual(values = c("#3FBC73FF", "#3FBC73FF", "#31B57BFF", "#31B57BFF", "#440154FF", "#440154FF", "#46337EFF", "#46337EFF", "#3F4889FF", "#3F4889FF", "#365C8DFF", "#365C8DFF", "#2E6E8EFF", "#2E6E8EFF")) +
    scale_x_discrete(labels = samples) +
    theme_bw() +
    theme(text = element_text(size = 16, color = "black"),
          axis.text.x = element_text(size = 16, angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(size = 16, color = "black"),
          axis.title.x = element_text(size = 16, color = "black"),
          axis.title.y = element_text(size = 16,color = "black"),
          legend.position = "none") +
    labs(title = "") +
    labs( x = "", y = "tpm") +
    guides(fill = FALSE) +
    scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
    geom_hline(yintercept = median(log2(all_combined_tpm.df$TC_tpm +1)), linetype="dashed", color = "indianred1", size=1, alpha = 0.8)


  pdf("Figures/shifting_tpm_boxplots.pdf", height = 4, width = 4)
    print(p)
  dev.off()

#----select PN6, oocyte, 2-cell, 4-cell and all tc
  all_combined_tpm_sel.df <- all_combined_tpm.df[all_combined_tpm.df$sample == "PN6" | all_combined_tpm.df$sample == "PN14" | all_combined_tpm.df$sample == "oocyte" | all_combined_tpm.df$sample == "S2_cell" | all_combined_tpm.df$sample == "S4_cell" | all_combined_tpm.df$sample == "PN6_all_TCs" | all_combined_tpm.df$sample == "PN14_all_TCs" | all_combined_tpm.df$sample == "oocyte_all_TCs" | all_combined_tpm.df$sample == "S2_cell_all_TCs" | all_combined_tpm.df$sample == "S4_cell_all_TCs", ]

# drop unused sequence levels
  all_combined_tpm_sel.df$sample <- droplevels(all_combined_tpm.df$sample)

# set levels for plotting
  all_combined_tpm_sel.df$sample <- factor(all_combined_tpm_sel.df$sample, levels = c("PN6", "PN6_all_TCs", "PN14", "PN14_all_TCs", "oocyte", "oocyte_all_TCs", "S2_cell", "S2_cell_all_TCs", "S4_cell", "S4_cell_all_TCs"))

  col <- col_scheme[c("P6", "P14", "MII", "2-cell", "4-cell")]
  samples <- levels(all_combined_tpm_sel.df$sample)

# ggplot 
  p <- ggplot(all_combined_tpm_sel.df, aes(y = log2(TC_tpm +1), x = sample , fill = sample, alpha = 0.6)) +
    geom_boxplot(outlier.shape = NA, varwidth = T, notch = F)  +
    scale_fill_manual(values = c("#440154FF", "#440154FF", "#46337EFF", "#46337EFF", "#3F4889FF", "#3F4889FF", "#365C8DFF", "#365C8DFF", "#2E6E8EFF", "#2E6E8EFF")) +
    scale_x_discrete(labels = samples) +
    theme_bw() +
    theme(text = element_text(size = 16, color = "black"),
          axis.text.x = element_text(size = 16, angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(size = 16, color = "black"),
          axis.title.x = element_text(size = 16, color = "black"),
          axis.title.y = element_text(size = 16,color = "black"),
          legend.position = "none") +
    labs(title = "") +
    labs( x = "", y = "tpm") +
    guides(fill = FALSE) +
    scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
    geom_hline(yintercept = median(log2(all_combined_tpm.df$TC_tpm +1)), linetype="dashed", color = "indianred1", size=1, alpha = 0.8)


  pdf("Figures/shifting_tpm_oocyte_boxplots.pdf", height = 4, width = 4)
    print(p)
  dev.off()
