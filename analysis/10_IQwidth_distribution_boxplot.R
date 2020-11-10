#/*==========================================================================#*/
#' ## Extended Data Figure 2B
#/*==========================================================================#*/
# distribution of IQ-widths - boxplots

library(dplyr)

# load TC object
  tc_PGCs_oocyte_embryo_anno.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/tc_PGCs_oocyte_embryo_anno_grl.RDS")

  names <- names(tc_PGCs_oocyte_embryo_anno.grl)

# extract IQwidths
  IQ_width.dfl <- lapply(names, function(x) {
    return(data.frame(IQ_width = tc_PGCs_oocyte_embryo_anno.grl[[x]]$interquantile_width,
                      sample = rep(x, times = length(tc_PGCs_oocyte_embryo_anno.grl[[x]]))))
  })

# flatten list into a df
  IQ_width.df <- do.call(what = rbind, IQ_width.dfl)
  
# plotting IQ-widths as boxplots plots - distribution is too wide for violin plots
  col <- col_scheme

# rename labels
  samples <- c("E9.5", "E10.5", "E11.5", "E12.5F", "E12.5M", "E13.5F", "E13.5M", "E14.5F", "E14.5M", "E16.5F", "E16.5M", "P6", "P14", "MII", "2-cell", "4-cell")

p <- ggplot(IQ_width.df, aes(y = IQ_width, x = sample , fill = sample, alpha = 0.6)) +
  geom_boxplot(outlier.shape = NA, varwidth = T, notch = F) +
  scale_fill_manual(values = c("firebrick1", "#FDE725FF",  "#E8E419FF",  "#D1E11CFF",  "#B9DE28FF", 
                               "#A2DA37FF",  "#8BD646FF",  "#75D054FF",  "#61CB5FFF",  "#4FC46AFF", "#3FBC73FF",  "#31B57BFF",  "#440154FF", "#46337EFF", 
                               "#3F4889FF",  "#365C8DFF",  "#2E6E8EFF")) +
  scale_x_discrete(labels = samples) +
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "none") +
  labs(title = "") +
  labs( x = "", y = "IQ-width") +
  guides(fill = FALSE) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, 20)) +
  geom_hline(yintercept = median(IQ_width_sel.df$IQ_width), linetype="dashed", color = "indianred1", size=1, alpha = 0.8)


pdf("Figures/boxplot_iq_width_pgc_oocyte.pdf", height = 5, width = 5)
  print(p)
dev.off() 