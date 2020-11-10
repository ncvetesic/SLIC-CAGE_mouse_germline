#/*==========================================================================#*/
#' ## Figure 2B
#/*==========================================================================#*/
# Distribution of distances between shifted dominant TSSs
# Shifting promoters have to be called first and domTSS distances calculated - see shifting_promoters.R

  domTSS_all_vs_E14_mESC_shift.grl <- readRDS("merged_replicates/intermediate_data/domTSS_all_vs_E14_mESC_shift_grl.RDS")

# plot distribution of domTSS switch 
  domTSS_dist.l <- lapply(domTSS_all_vs_E14_mESC_shift.grl, function(x) {
    df <- data.frame(domTSS_dist = as.numeric(x$domTSS_dist))
  })

# add sample name
  for(i in 1:length(domTSS_dist.l)) {
    domTSS_dist.l[[i]]$sample <- rep(names(domTSS_dist.l)[i], times = nrow(domTSS_dist.l[[i]]))
  }

# flatten list to a dataframe
  domTSS_dist.df <- do.call(rbind, domTSS_dist.l)

# change levels
  domTSS_dist.df$sample <- factor(domTSS_dist.df$sample, levels = names(domTSS_all_vs_E14_mESC_shift.grl))

# plot distribution of domTSS distances 
  col = rep("dodgerblue", times = 14)

  p <- ggplot(domTSS_dist.df, aes(x = domTSS_dist,  fill = sample), alpha = 0.5) +
    geom_histogram (binwidth = 1, col = "black", lwd = 0.125) +
    scale_fill_manual("sample", values = col) +
    theme_bw() +
    theme(text = element_text(size = 14, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(y = "Frequency", x = NULL) + 
    coord_cartesian( xlim = c(-100, 100), ylim = c(0, 250)) +
    guides(fill = FALSE)

  pdf("merged_replicates/results/shifting_promoters/mESC/distribution_shift_domTSS.pdf", width = 12, height = 12)
    p + facet_wrap(~ sample, ncol = 3)
  dev.off()


# print zoom in
  p <- ggplot(domTSS_dist.df, aes(x = domTSS_dist,  fill = sample), alpha = 0.5) +
    geom_histogram(binwidth = 1, col = "black", lwd = 0.125) +
    scale_fill_manual("sample", values = col) +
    theme_bw() +
    theme(text = element_text(size = 14, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(y = "Frequency", x = NULL) + 
    coord_cartesian( xlim = c(-100, 100), ylim = c(0, 50)) +
    guides(fill = FALSE)

  pdf("merged_replicates/results/shifting_promoters/mESC/distribution_shift_domTSS_zoom.pdf", width = 12, height = 12)
    p + facet_wrap(~ sample, ncol = 3)
  dev.off()