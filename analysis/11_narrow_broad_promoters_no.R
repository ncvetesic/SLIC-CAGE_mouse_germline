#/*==========================================================================#*/
#' ## Extended Data Figure 2C
#/*==========================================================================#*/
# number of narrow and broad promoters identified in each sample

# import TC data
  tc_PGCs_oocyte_embryo_anno.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/tc_PGCs_oocyte_embryo_anno_grl.RDS")
  
# exclude weird chromosomes
  tc_PGCs_oocyte_embryo_anno.grl <- lapply(tc_PGCs_oocyte_embryo_anno.grl, function(x){
    gr <- dropSeqlevels(x, "chrM", pruning.mode = "coarse")
  })

# divide into narrow/sharp and broad - add mcols column
  tc_PGCs_oocyte_embryo_anno_broadness.grl <- lapply(tc_PGCs_oocyte_embryo_anno.grl, function(x) {
    broadness <- ifelse(x$interquantile_width < 9, "sharp", "broad")
    x$broadness <- broadness
    return(x)
  })

# broadness stats
  broadness_stats.dfl <- lapply(1:length(tc_PGCs_oocyte_embryo_anno_broadness.grl), function(x) {
    df <- data.frame(table(tc_PGCs_oocyte_embryo_anno_broadness.grl[[x]]$broadness))
    colnames(df) <- c("broadness", "number")
    df$sample <- names(tc_PGCs_oocyte_embryo_anno_broadness.grl)[x]
    return(df)
  })

# add names to elements
  names(broadness_stats.dfl) <- names(tc_PGCs_oocyte_embryo_anno.grl)

# flatten to a single dataframe
  broadness_stats.df <- do.call(rbind, broadness_stats.dfl)

# change levels
  broadness_stats.df$sample <- droplevels(broadness_stats.df$sample)
  broadness_stats.df$broadness <- factor(broadness_stats.df$broadness, levels = c("sharp", "broad"))
  broadness_stats.df$sample <- factor(broadness_stats.df$sample, levels = names(tc_PGCs_oocyte_embryo_anno_broadness.grl))

# ggplot plotting
  library(ggplot2)
  p <- ggplot(broadness_stats.df, aes(x = sample, y = number, fill = broadness), alpha = 0.7) +
    geom_bar(stat = "identity", width = 0.75, colour = "black", lwd = 0.5) +
    scale_fill_manual("Features", values = c("darkgoldenrod1", "dimgray")) +
    theme_light() +
    theme(text = element_text(size = 16, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black",),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(colour = "black")) +
    labs(y = "Number", x = NULL) +
    scale_y_continuous(limits = c(0, 40000), breaks = seq(0, 40000, by = 5000))

  p <- p + geom_text(aes(label = number), size = 3, hjust = 0.5, vjust = 3, position = "stack")

# print in pdf
  pdf("Figures/sharp_broad_barplot.pdf", width = 10, height = 6)
    print(p)
  dev.off()