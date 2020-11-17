#/*==========================================================================#*/
#' ## Extended Data Figure 11A,B
#/*==========================================================================#*/
# Somatic CTSSs expression levels

# need SOM clustering of dominant CTSSs first - see Extended Data Figure 15a

# using som clusters of dominant CTSSs 17 and 18 - somatic
  library(ggplot2)

# load SOM divided dominant CTSS object
  domCTSSs_SOM25_PGC_F_oocyte_promOnly.grl <- readRDS("intermediate_data/domCTSSs_SOM25_PGC_F_oocyte_promOnly.grl.RDS")

# add SOM class column to each GRanges list element  
  CTSSs_SOMclass.grl <- GRangesList(domCTSSs_SOM25_PGC_F_oocyte_promOnly.grl)
  CTSSs_SOMclass.grl <- lapply(1:length(CTSSs_SOMclass.grl), function(x) {
    CTSSs_SOMclass.grl[[x]]$SOM_class <- names(CTSSs_SOMclass.grl)[x]
    return(CTSSs_SOMclass.grl[[x]])
  })

# unlist or concatenate to one GRanges
  CTSSs_SOMclass.gr <- do.call(c, unlist(CTSSs_SOMclass.grl))

# separate into grangeslist by samples
  samples <- colnames(mcols(CTSSs_SOMclass.gr))[-c(1:3, 11, 12)]

  CTSSs_SOMclass_samples.grl <- lapply(1:length(samples), function(x) {
    gr <- CTSSs_SOMclass.gr
    mcols(gr) <- NULL
    gr$tpm <- mcols(CTSSs_SOMclass.gr)[, samples[x]]
    gr$SOM_class <- mcols(CTSSs_SOMclass.gr)$SOM_class
    gr$CTSSid  <- mcols(CTSSs_SOMclass.gr)$CTSSid 
    return(gr)
  })

# add sample names
  names(CTSSs_SOMclass_samples.grl) <- samples

# extract som 17 and 18 from E16.5F, P6, P14 and MII plot violin plot
  selected_SOM17_18_ctss.grl <- CTSSs_SOMclass_samples.grl[-c(1:3)]

# extract 17 and 18
  selected_SOM17_18_ctss.grl <- lapply(selected_SOM17_18_ctss.grl, function(x) {
    gr <- x[mcols(x)$SOM_class == "som17" | mcols(x)$SOM_class == "som18"]
    return(gr)
  })

# add sample name to each
  selected_SOM17_18_ctss.grl <- lapply(1:length(selected_SOM17_18_ctss.grl), function(x) {
    selected_SOM17_18_ctss.grl[[x]]$sample <- names(selected_SOM17_18_ctss.grl)[x]
    return(selected_SOM17_18_ctss.grl[[x]])
  })

# flatten to dataframe and unlist
  selected_SOM17_18_ctss.dfl <- lapply(selected_SOM17_18_ctss.grl, function(x) {
    df <- as.data.frame(x)
    return(df)
  })

  selected_SOM17_18_ctss.df <- do.call(rbind, selected_SOM17_18_ctss.dfl)

# selected columns
  selected_SOM17_18_ctss_sel.df <- selected_SOM17_18_ctss.df[, 6:9]

# convert to factors
  selected_SOM17_18_ctss_sel.df$SOM_class <- factor(selected_SOM17_18_ctss_sel.df$SOM_class, levels = c("som17", "som18"))
  selected_SOM17_18_ctss_sel.df$sample <- factor(selected_SOM17_18_ctss_sel.df$sample, levels = c("E16_5F", "PN6", "PN14", "oocyte"))

# plotting
  library(ggplot2)
  col <- col_scheme[c("E16.5F", "P6", "P14", "MII")]
  names(col) <- c("E16_5F", "PN6", "PN14", "oocyte")
  med_value <- median(log2(selected_SOM17_18_ctss_sel.df$tpm +1))

  p <-  ggplot(selected_SOM17_18_ctss_sel.df, aes(x = sample, y = log2(tpm +1), fill = sample), inherit.aes = FALSE) + 
    geom_violin(alpha = .4, trim = FALSE, lwd = 0.25) +
    geom_boxplot(width=0.1, alpha= 0.6, outlier.alpha = 0.2) +
    scale_fill_manual(values = col) +
    theme_bw() +
    theme(text = element_text(size = 16),
          axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 1, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          strip.background =element_rect(fill="white")) +
    coord_cartesian(ylim=(c(0, 6))) +
    geom_hline(yintercept = med_value, linetype="dashed", color = "indianred1", size=1, alpha = 0.8)
  
  
  pdf("Figures/SOM_domCTSSs_somatic_in_oocyte.pdf", height = 4, width = 8)
    p + facet_wrap(~ SOM_class)
  dev.off() 

# extract same CTSSs in P14 and P6 and show tpm scatterplot
  selected_SOM17_18_ctss_P6_P14.df <- selected_SOM17_18_ctss_sel.df[selected_SOM17_18_ctss_sel.df$sample == "PN6" | selected_SOM17_18_ctss_sel.df$sample == "PN14", ]

# drop levels
  selected_SOM17_18_ctss_P6_P14.df$sample <- droplevels(selected_SOM17_18_ctss_P6_P14.df$sample)

# separate into 2 dataframes
  P6.df <- selected_SOM17_18_ctss_P6_P14.df[selected_SOM17_18_ctss_P6_P14.df$sample == "PN6", ]
  
  P14.df <- selected_SOM17_18_ctss_P6_P14.df[selected_SOM17_18_ctss_P6_P14.df$sample == "PN14", ]

# merge back datasets by common CTSSid
  #rename columns
    names(P6.df) <- paste0(colnames(P6.df), "_P6")
    names(P14.df) <- paste0(colnames(P14.df), "_P14")

    P6_P14_merge.df <- merge(P6.df, P14.df, by.x = "CTSSid_P6", by.y = "CTSSid_P14")


# plotting  
  p <- ggplot(P6_P14_merge.df, aes(x = log2(P6_P14_merge.df$tpm_P6 + 1), y = log2(P6_P14_merge.df$tpm_P14 + 1))) +
    geom_point(shape = 21, size = 2, alpha = 0.15, shape = 21, fill = "black") +
    theme_light() +
    theme(text = element_text(size = 16, colour = "black"),
          axis.text.x = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title.x = element_text(size = 16, colour = "black"),
          axis.title.y = element_text(size = 16, colour = "black"),
          strip.background =element_rect(fill="white")) +
    geom_abline(x = 0:20, y = 0:20, lty = "dashed", colour = "red", size = 0.8) +
    scale_y_continuous(limits = c(0, 3)) +
    scale_x_continuous(limits = c(0, 3)) +
    geom_density_2d(colour = "darkorange") +
    coord_fixed(ratio = 1)

  pdf("Figures/SOM_domCTSSs_somatic_tpm_scatter.pdf", height = 4, width = 8)
    p + facet_wrap(~ SOM_class_P6)
  dev.off()