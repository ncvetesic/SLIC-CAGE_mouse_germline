#/*==========================================================================#*/
#' ## Figure 1C
#/*==========================================================================#*/

library(corrplot)

# extract TPM values normalized for all samples
  tag.count <- CAGEset_PGC_embryo_merged@normalizedTpmMatrix

#calculate correlations
  corr_mat <- cor(as.matrix(tag.count), method = "pearson")
  names <- colnames(tag.count)

  rownames(corr_mat) <- names
  colnames(corr_mat) <- names

# ploting final graph
  library(RColorBrewer)
  col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                             "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                             "#4393C3", "#2166AC", "#053061"))

  quartz(type = "pdf", file = "Figures/corrPlot_selected_merged_full_nonum.pdf", width = 10 , height = 10, family = "Helvetica", dpi = 300)
  corrplot(corr_mat, 
           method = "color", 
           tl.pos = "td",
           col = col2(20),
           cl.lim = c(0, 1),
           tl.col = "black",
           tl.srt = 45,
           tl.cex = 1,
           outline = T, 
           addgrid.col = "darkgray",
           bg = "white",
           type = "upper",
           is.corr = F)
  dev.off()