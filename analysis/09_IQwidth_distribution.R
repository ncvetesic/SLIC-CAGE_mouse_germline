#/*==========================================================================#*/
#' ## Extended Data Figure 1A
#/*==========================================================================#*/
# distribution of IQ-widths

# extract sample names
  all_samples <- sampleLabels(CAGEset_PGC_embryo_merged)

# extract interquantile widths
  iq.l <- lapply(tc_PGCs_oocyte_embryo_anno.grl, function(x) data.frame(iq_width = x$interquantile_width))

# add names to as a column to each sample
  for (i in 1:length(all_samples)) {
    iq.l[[i]]$sample <- rep(all_samples[[i]], nrow(iq.l[[i]]))
  }

# combine to 1 dataframe
  iq.df <- do.call(rbind, iq.l)

# set levels for plotting - sample levels
  iq.df$sample <- factor(iq.df$sample, levels = all_samples)

# plotting
  library(ggplot2)
  p <- ggplot(iq.df, aes(iq_width)) +
    geom_histogram(aes(x = iq_width, y = ..ncount..* 100), binwidth = 1, 
                   fill = "gray60", col = "black", size = 0.1) +
    theme_bw() +
    theme(text = element_text(size = 16, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black", size = 14),
          axis.text.y = element_text(colour = "black", size = 14),
          strip.background = element_rect(colour = "black", fill = "gray87")) +
    labs(x = "Interquantile width", y = "Percentage") +
    coord_equal(ratio = 1) +
    xlim(0, 150)

pdf(file = "Figures/supplementary/IQ_width_all.pdf", height = 8, width = 12)
  p + facet_wrap(~ sample, ncol = 6)
dev.off()
