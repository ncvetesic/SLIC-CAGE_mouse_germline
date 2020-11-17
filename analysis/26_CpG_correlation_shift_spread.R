#/*==========================================================================#*/
#' ## Extended Data Figure 12H
#/*==========================================================================#*/
# Correlation of tag cluster signal spread width in E16.5F vs mESC shifting promoters and its CpG content

### Spread vs CpG E16_5F diff broad and narrow separate
  # broad in E16.5F regions
  # select regions of spreading - both ranges have to be converted into a GRangesList as each range difference can produce more regions of difference - hence it has to be a list for each pair we are comparing
    mESC_vs_E16_5F_broader.gr <- as(mESC_vs_E16_5F_broader.gr, "GRangesList")
    E16_5F_TCs_mESC_vs_E16_5F_broader.gr <- as(E16_5F_TCs_mESC_vs_E16_5F_broader.gr, "GRangesList")
    
    E16_5_vs_mESC_setdiff.grl <- setdiff(E16_5F_TCs_mESC_vs_E16_5F_broader.gr, mESC_vs_E16_5F_broader.gr)

  # unlist
    E16_5_vs_mESC_setdiff.gr <- unlist(E16_5_vs_mESC_setdiff.grl)

  # width
  E16_5_vs_mESC_setdiff_width <- width(E16_5_vs_mESC_setdiff.gr)

  # extract sequence
    E16_5_vs_mESC_setdiff_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, E16_5_vs_mESC_setdiff.gr)

  # count CpG in each
    E16_5_vs_mESC_setdiff_count <- dinucleotideFrequency(E16_5_vs_mESC_setdiff_seq, step = 1, as.prob = F)
    E16_5_vs_mESC_setdiff_CpG_count <- E16_5_vs_mESC_setdiff_count[, "CG"] + E16_5_vs_mESC_setdiff_count[, "GC"]

  # add as metacolumn
    E16_5_vs_mESC_setdiff.gr$CpG_count <- E16_5_vs_mESC_setdiff_CpG_count

  # normalise by length
    E16_5_vs_mESC_setdiff_CpG_count_norm <- E16_5_vs_mESC_setdiff_CpG_count/(E16_5_vs_mESC_setdiff_width-1)

  # plot correlation of width and CpG count
    corr_CpG_width_diff.df <- data.frame("width" = E16_5_vs_mESC_setdiff_width, "CpG_norm_count" = E16_5_vs_mESC_setdiff_CpG_count_norm)

  # remove NA values - in normalised CpG count
    corr_CpG_width_diff.df <- corr_CpG_width_diff.df[!is.na(corr_CpG_width_diff.df$CpG_norm_count), ]
    corr_CpG_width_diff.df$width <- as.numeric(corr_CpG_width_diff.df$width)

  # add spearman correlation
    p <- ggplot(corr_CpG_width_diff.df) +
      geom_point(shape = 21, aes(x = CpG_norm_count, y = width), size = 3, alpha = 0.6, shape = 21, fill = "black") +
      theme_light() +
      theme(text = element_text(size = 20, colour = "black"), 
            axis.text.x = element_text(size = 18, colour = "black"),
            axis.text.y = element_text(size = 18, colour = "black"),
            legend.title = element_blank()) +
      guides(colour = guide_legend(reverse=T)) +
      scale_y_continuous(name = "Spread width", limits = c(0, 200), seq(0, 200, 25)) +
      scale_x_continuous(name = "CpG count norm", limits = c(0, 0.8), seq(0, 0.8, 0.1)) +
      geom_smooth(method = 'loess', span = 0.75, aes(x = CpG_norm_count, y = width), lty = "dashed", colour = "red", size = 1)
    
      pdf("Figures/spread_cpg_corr_braoder_loess.pdf", width = 6, height = 6)
        print(p)
      dev.off() 

# narrower in E16.5F regions
# select regions of spreading - both ranges have to be converted into a GRangesList as each range difference can produce more regions of difference - hence it has to be a list for each pair we are comparing
  mESC_vs_E16_5F_narrower.gr <- as(mESC_vs_E16_5F_narrower.gr, "GRangesList")
  E16_5F_TCs_mESC_vs_E16_5F_narrower.gr <- as(E16_5F_TCs_mESC_vs_E16_5F_narrower.gr, "GRangesList")
  
  E16_5_vs_mESC_setdiff.grl <- setdiff(E16_5F_TCs_mESC_vs_E16_5F_narrower.gr, mESC_vs_E16_5F_narrower.gr)

# unlist
  E16_5_vs_mESC_setdiff.gr <- unlist(E16_5_vs_mESC_setdiff.grl)

# width
  E16_5_vs_mESC_setdiff_width <- width(E16_5_vs_mESC_setdiff.gr)

# extract sequence
  E16_5_vs_mESC_setdiff_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, E16_5_vs_mESC_setdiff.gr)

# count CpG in each
  E16_5_vs_mESC_setdiff_count <- dinucleotideFrequency(E16_5_vs_mESC_setdiff_seq, step = 1, as.prob = F)
  E16_5_vs_mESC_setdiff_CpG_count <- E16_5_vs_mESC_setdiff_count[, "CG"] + E16_5_vs_mESC_setdiff_count[, "GC"]

# add as metacolumn
  E16_5_vs_mESC_setdiff.gr$CpG_count <- E16_5_vs_mESC_setdiff_CpG_count

# normalise by length
  E16_5_vs_mESC_setdiff_CpG_count_norm <- E16_5_vs_mESC_setdiff_CpG_count/(E16_5_vs_mESC_setdiff_width-1)

# plot correlation of width and CpG count
  corr_CpG_width_diff.df <- data.frame("width" = E16_5_vs_mESC_setdiff_width, "CpG_norm_count" = E16_5_vs_mESC_setdiff_CpG_count_norm)

# remove NA values - in normalised CpG count
  corr_CpG_width_diff.df <- corr_CpG_width_diff.df[!is.na(corr_CpG_width_diff.df$CpG_norm_count), ]
  corr_CpG_width_diff.df$width <- as.numeric(corr_CpG_width_diff.df$width)

# spearman correlation
  p <- ggplot(corr_CpG_width_diff.df) +
    geom_point(shape = 21, aes(x = CpG_norm_count, y = width), size = 3, alpha = 0.6, shape = 21, fill = "black") +
    theme_light() +
    theme(text = element_text(size = 20, colour = "black"), 
          axis.text.x = element_text(size = 18, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black"),
          legend.title = element_blank()) +
    guides(colour = guide_legend(reverse=T)) +
    scale_y_continuous(name = "Spread width", limits = c(0, 200), seq(0, 200, 25)) +
    scale_x_continuous(name = "CpG count norm", limits = c(0, 0.8), seq(0, 0.8, 0.1)) +
    geom_smooth(method = 'loess', span = 0.75, aes(x = CpG_norm_count, y = width), lty = "dashed", colour = "red", size = 1)
  
  pdf("Figures/spread_cpg_corr_narrower_loess.pdf", width = 6, height = 6)
    print(p)
  dev.off()  

# plot using ggpubr
  library(ggpubr)
  
  pdf("Figures/spread_cpg_corr_narrower.pdf", width = 6, height = 6)
  ggscatter(corr_CpG_width_diff.df, x = "CpG_norm_count", y = "width", add = "reg.line",
            conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman", cor.coef.coord	= c(0.0, 200),
            
            xlab = "CpG count norm", ylab = "Spread width", size = 3, alpha = 0.6, shape = 21, 
            fill = "lightgray", xlim = c(0, 0.4), ylim = c(0, 200),
            color = "black", ggtheme = theme_light(),
            add.params = list(color = "red", lty = "dashed", size = 1))
  dev.off()