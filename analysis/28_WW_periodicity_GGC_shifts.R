#/*==========================================================================#*/
#' ## Extended Data Figure 12J, K
#/*==========================================================================#*/
# WW periodicity in E16.5F vs mESC shifting promoters
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(seqPattern)

## Metaplot data function from heatmap object
   metaplot_data <- function(hm.l, motif = FALSE, up = 500,  down = 500, smooth_span = 3, smooth_step = 1) {
    
    # path - sets where to save the image and image name
    # hm.l is a heatmap object (not a list) - Malcolm's heatmap object with extractable image
    # classes is a factor  class that divides a metaplot
    
    # subset matrix in classes
    matrix_l <- image(hm.l)
    
    # metaplot calculation and smoothing
    library(zoo)
    matrix_sum.l <-  colSums(matrix_l)/nrow(matrix_l)
    df_smooth.l <- rollapply(data = matrix_sum.l, FUN = mean, width = smooth_span, by = smooth_step)
    
    df_smooth.df <- as.data.frame(df_smooth.l)
    
    # add coordinate column
    if (motif == FALSE) {
      df_smooth.df$coord <- rollapply(data = -up:(down - 1), FUN = mean, width = smooth_span, by = smooth_step)
    } else {
      df_smooth.df$coord <- rollapply(data = -up:(down - (up + down - nrow(df_smooth.df) - 1)), FUN = mean, width = smooth_span, by = smooth_step)
    }
    
    # format dataframe to ggplot acceptable
    library(tidyr)
    #df_smooth.gg <- gather(df_smooth.df, 
    #... = colnames(df_smooth.df[1:length(levels(classes))]),
    #key = "type",
    #value = "value")
    #df_smooth.gg$type <- factor(df_smooth.gg$type, levels = levels(classes))
    
    return(df_smooth.df)
  }
  
# load domTSS shifting tag clusters
  domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/domTSS_PGC_oocyte_embryo_vs_mESC_shift_grl.RDS")

# centred on X are centred on NON-mESC samples
  domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/domTSS_PGC_oocyte_embryo_vs_mESC_shift_shift_X_grl.RDS")

  domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl <- lapply(domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl, function(x){
    gr <- dropSeqlevels(x, "chrM", pruning.mode = "coarse")
  })

  domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl <- lapply(domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl, function(x){
    gr <- dropSeqlevels(x, "chrM", pruning.mode = "coarse")
  })

  samples <- names(domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl)

# attach iq_width to switching promoters
  tc_PGCs_oocyte_embryo_anno.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/tc_PGCs_oocyte_embryo_anno_grl.RDS")

# exclude mESC
  tc_PGCs_oocyte_embryo_anno_nomESC.grl <- tc_PGCs_oocyte_embryo_anno.grl[-1]

# overlap shifting promoters and tc to get iq-width/annotate with mcols
  domTSS_overlap.l <- list()
  for(i in 1:length(domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl)) {
    overlap <- findOverlapPairs(domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl[[i]], tc_PGCs_oocyte_embryo_anno_nomESC.grl[[i]])
    domTSS_overlap.l[i] <- overlap@first
    mcols(domTSS_overlap.l[[i]]) <- cbind(mcols(overlap@first), mcols(overlap@second))
  }
  names(domTSS_overlap.l) <- names(domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl)

# filter to include shifting promoters that moved by at least one bp
  domTSS_overlap.l <- lapply(domTSS_overlap.l, function(x) {
    sample <- x[abs(x$domTSS_dist) >= 1]
    return(sample)
  })

# create windows centered on domTSS -100, +300
  up <- 100
  down <- 300
  range <- c(-up, down)
  win <- up + down

  domTSS_E14_win.grl <- sapply(domTSS_overlap.l, function(x) promoters(x, up = up, down = down))

# remove out of bound ranges
  domTSS_E14_win.grl <- sapply(domTSS_E14_win.grl, function(x) x[width(trim(x)) == win])

# extract sequence 
  domTSS_E14_seq.l <- sapply(domTSS_E14_win.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))

# use malcolms heatmaps to create smoothed pattern heatmaps and sum it into a metaplot..- plot using ggplot2
  library(heatmaps)

# plot dinucleotide plots centered on dominant TSS
  pattern <- "WW"
  hm.l <- list()
  hm_smoothed.l <- list()

# sort index - by distance
  sort.l <- lapply(domTSS_E14_win.grl, function(x) {
    return(order(x$domTSS_dist, decreasing = T))
  })

# plot metaplots 
# produce heatmaps prior to plotting (no smoothing)
  hm_broad.l <- list()

  for(i in 1:length(domTSS_E14_seq.l)) {
    hm_broad.l[[i]] <- PatternHeatmap(domTSS_E14_seq.l[[i]], pattern = pattern, coords = range, label = paste(pattern, ", ", samples[i], sep = ""))
  }
  names(hm_broad.l) <- names(domTSS_E14_seq.l)

# convert heatmaps to metaplot data - returns df 
  meta_broad.l <- lapply(hm_broad.l, function(x) metaplot_data(x, up = 100, down = 300))

# create a list of dataframes (with sharp and broad column)
  all_df.l <- list()

  samples <- names(domTSS_E14_seq.l)

  for (i in 1:length(samples)) {
    all_df.l[[i]] <- data.frame("x_coord" = meta_broad.l[[i]]$coord, 
                                "broad_occurrence" = meta_broad.l[[i]]$df_smooth.l,
                                "sample" = rep(samples[i], times = length(meta_broad.l[[1]]$df_smooth.l)))
  }
  names(all_df.l) <- samples

# remove tbp2 KO
  all_df.l <- all_df.l[-c(14, 16)]

# collapse the list into a dataframe
  all_df <- do.call(rbind, all_df.l)
  rownames(all_df) <- 1:nrow(all_df)

# plot metaplots using ggplot2
library(ggplot2)
  p <- ggplot(all_df) +
    geom_line(lty = 1, aes(x = x_coord, y = broad_occurrence, colour = sample), size = 0.6) +
    scale_color_manual(values = c("#FDE725FF", "#E8E419FF", "#D1E11CFF", "#B9DE28FF", "#A2DA37FF",
                                  "#8BD646FF", "#75D054FF", "#61CB5FFF", "#4FC46AFF", "#3FBC73FF", "#31B57BFF",
                                  "#440154FF", "#481B6DFF",  "#46337EFF",  "#3F4889FF",  "#365C8DFF",  "#2E6E8EFF"), 
                       name = NULL) +
    theme_light()   +
    theme(text = element_text(size = 20, colour = "black"),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          legend.position = "right",
          legend.title = element_blank()) +
    guides(colour = guide_legend(reverse=F)) +
    geom_vline(xintercept = 0, lty = "dashed", colour = "gray48", size = 0.6) +
    scale_y_continuous(name = "Relative frequency", limits = c(0, 0.5), breaks = c(0, 0.25, 0.5)) +
    scale_x_continuous(name = "Distance to dominant TSS (bp)", limits = c(-100, 300), seq(-100, 300, by = 100))
  
  p <- p + facet_wrap(~ sample, ncol = 3)

  pdf("Figures/WW_domTSS_dist1.pdf", height = 12, width = 11)
    print(p)
  dev.off()

# plot selected samples
  all_df_sub.gg <- all_df[all_df$sample %in% c("E16_5F", "PN6"), ]
  all_df_sub.gg$sample <- droplevels(all_df_sub.gg$sample)

  library(ggplot2)
  p <- ggplot(all_df_sub.gg) +
    geom_line(lty = 1, aes(x = x_coord, y = broad_occurrence, colour = sample), size = 0.6) +
    scale_color_manual(values = c("#3FBC73FF", "#440154FF"), name = NULL) +
    theme_light()   +
    theme(text = element_text(size = 20, colour = "black"),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          legend.position = "right",
          legend.title = element_blank()) +
    guides(colour = guide_legend(reverse=F)) +
    geom_vline(xintercept = 0, lty = "dashed", colour = "gray48", size = 0.6) +
    scale_y_continuous(name = "Relative frequency", limits = c(0, 0.5), breaks = c(0, 0.25, 0.5)) +
    scale_x_continuous(name = "Distance to dominant TSS (bp)", limits = c(-100, 300), seq(-100, 300, by = 100))
  
  pdf("Figures/WW_domTSS_dist_sel_split.pdf", height = 4, width = 11)
    print(p + facet_wrap(~ sample, ncol = 2))
  dev.off()
  
  pdf("Figures/WW_domTSS_dist_sel.pdf", height = 4, width = 6)
    print(p)
  dev.off()

# zoomed metaplots for the inset = only broad pattern occurence
  library(ggplot2)
  p <- ggplot(all_df, aes(x = x_coord, y = broad_occurrence, colour = sample)) +
    geom_line(lty = 1, size = 0.6, legend = FALSE) +
    scale_color_manual(values = c("#FDE725FF", "#E8E419FF", "#D1E11CFF", "#B9DE28FF", "#A2DA37FF",
                                  "#8BD646FF", "#75D054FF", "#61CB5FFF", "#4FC46AFF", "#3FBC73FF", "#31B57BFF",
                                  "#440154FF", "#481B6DFF",  "#46337EFF",  "#3F4889FF",  "#365C8DFF",  "#2E6E8EFF"), 
                       name = NULL)  +
    theme_light()   +
    theme(text = element_text(size = 20, colour = "black"),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          legend.position = "right",
          legend.title = element_blank()) +
    guides(colour = guide_legend(reverse=F)) +
    coord_cartesian(xlim = c(50, 200), ylim = c(0.05, 0.2))
  p <- p + facet_wrap(~ sample, ncol = 4)
  
  pdf("Figures/WW_domTSS_dist_zoom.pdf", height = 12, width = 11)
    print(p)
  dev.off()

# plot zoomed selected samples
# plot selected samples
  all_df_sub.gg <- all_df[all_df$sample %in% c("E16_5F", "PN6"), ]
  all_df_sub.gg$sample <- droplevels(all_df_sub.gg$sample)

  library(ggplot2)
  p <- ggplot(all_df_sub.gg) +
    geom_line(lty = 1, aes(x = x_coord, y = broad_occurrence, colour = sample), size = 0.6) +
    scale_color_manual(values = c("#3FBC73FF", "#440154FF"), name = NULL) +
    theme_light()   +
    theme(text = element_text(size = 20, colour = "black"),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          legend.position = "right",
          legend.title = element_blank()) +
    guides(colour = guide_legend(reverse=F)) +
    geom_vline(xintercept = 0, lty = "dashed", colour = "gray48", size = 0.6) +
    scale_y_continuous(name = "Relative frequency") +
    scale_x_continuous(name = "Distance to dominant TSS (bp)") +
    coord_cartesian(xlim = c(50, 200), ylim = c(0.075, 0.15))

  pdf("Figures/WW_domTSS_dist_sel_split_zoom.pdf", height = 4, width = 11)
    print(p + facet_wrap(~ sample, ncol = 2))
  dev.off()

  pdf("Figures/WW_domTSS_dist_sel_zoom.pdf", height = 4, width = 6)
    print(p)
  dev.off()

#------ check nucleosome positioning signal in windows centered on randomly selected true CTSS ------#
#extract CTSS dataframe with normalized TPM per CTSS per sample
  CTSSnorm <- CTSSnormalizedTpm(CAGEset_PGC_embryo_merged)

# filter to have at least 1TPM 
  filtr_idx.list <- list()
  samples_all <- colnames(CTSSnorm)[-c(1:3)]
  
  for (i in 1:length(samples_all)) {
    filtr_idx.list[[i]] <- c(CTSSnorm[, samples_all[i]] >= 1)
  }
  names(filtr_idx.list) <- samples_all

# create GRanges object CTSSs (filter per sample to have >= 1 TPM)
  filtr.ctss.grl <- lapply(filtr_idx.list, function(x) GRanges(seqnames = CTSSnorm[x, ]$chr, 
                                                               ranges = IRanges(start = CTSSnorm[x, ]$pos, width = 1),
                                                               strand = CTSSnorm[x, ]$strand))

# create granges object of shifting promoters filtered at tc level
  TC_overlap.l <- lapply(domTSS_overlap.l, function(x) {
    gr <- GRanges(seqnames = seqnames(x),
                  ranges = IRanges(start = x$q_0.1,
                                   end = x$q_0.9),
                  strand = strand(x))
    mcols(gr) <- mcols(x)
    return(gr)
  })

# overlap tag clusters and CTSS - to connect which CTSSs are in which tagClusters, but only CTSSs with >= 1 TPM (not all CTSS)
  tc_CTSS_idx.l <- list()

  for (i in 1:length(TC_overlap.l)) {
    tc_CTSS_idx.l[[i]] <- findOverlaps(TC_overlap.l[[i]], filtr.ctss.grl[[i]], select = "all")
    # convert to a dataframe
    tc_CTSS_idx.l[[i]] <- as.data.frame(tc_CTSS_idx.l[[i]])
  }
  names(tc_CTSS_idx.l) <- names(TC_overlap.l)

# for each tag cluster select a random real CTSS to use as dominant - end object should be of equal length as each tag cluster
# collapse by tag cluster number
# create an index of tag clusters that have overlapping CTSSs with tpm >=1 
  tc.idx.l <- lapply(tc_CTSS_idx.l, function(x) unique(x$queryHits))

# sample random CTSS per overlapping tag cluster 
# select random CTSS per tag cluster - take indexes of CTSSs that overlap tag clusters, and sample one of them
  set.seed(10)
  random.ctss.l <- list()
  for (i in 1:length(tc.idx.l)) {
    random.ctss.l[[i]] <- sapply(tc.idx.l[[i]], function(x) sample(tc_CTSS_idx.l[[i]][tc_CTSS_idx.l[[i]]$queryHits == x, "subjectHits"], 1))
  }
  names(random.ctss.l) <- names(tc_CTSS_idx.l)

# subset tag-clusters in each sample based on the ones that where overlapped with ctss's
  tc_ctss_sub.l <- list()

  for (i in 1:length(tc.idx.l)) {
    tc_ctss_sub.l[[i]] <- TC_overlap.l[[i]][tc.idx.l[[i]]]
  }
  names(tc_ctss_sub.l) <- names(TC_overlap.l)

# get random CTSS coordinates (I have indices)
  random.ctss.coord.l <- list()

  for(i in 1:length(random.ctss.l)) {
    random.ctss.coord.l[[i]] <- filtr.ctss.grl[[i]][random.ctss.l[[i]], ]
  }

# attach CTSS information to corresponding tag cluster
  for (i in 1:length(tc_ctss_sub.l)) {
    tc_ctss_sub.l[[i]]$random_ctss_coord <- start(random.ctss.coord.l[[i]])
    tc_ctss_sub.l[[i]]$random_ctss_chr <- chrom(random.ctss.coord.l[[i]])
    seqinfo(random.ctss.coord.l[[i]]) <- seqinfo(tc_ctss_sub.l[[i]])[seqlevels(random.ctss.coord.l[[i]])]
    random.ctss.coord.l[[i]]$interquantile_width <- tc_ctss_sub.l[[i]]$interquantile_width
  }
  names(random.ctss.coord.l) <- names(tc_ctss_sub.l)

# # create windows centered on random TSS from that tag cluster -100, +300
  up <- 100
  down <- 300
  range <- c(-up, down)
  win <- up + down
  
  randomTSS_E14_win.grl <- sapply(random.ctss.coord.l, function(x) promoters(x, up = up, down = down))

# remove out of bound ranges
  randomTSS_E14_win.grl <- sapply(randomTSS_E14_win.grl, function(x) x[width(trim(x)) == win])

# check the length of the each range in the list
  sapply(randomTSS_E14_win.grl, length)

# extract sequence 
  randomTSS_E14_seq.l <- sapply(randomTSS_E14_win.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))

# exclude TBP2KO
  randomTSS_E14_seq.l <- randomTSS_E14_seq.l[-c(14, 16)]

# use malcolms heatmaps to create smoothed pattern heatmaps and sum it into a metaplot..- plot using ggplot2
  library(heatmaps)

# plot dinucleotide plots centered on dominant TSS
  pattern <- "WW"
  hm.l <- list()
  hm_smoothed.l <- list()

# plot metaplots 
# produce heatmaps prior to plotting (no smoothing)
  hm_broad.l <- list()

  for(i in 1:length(randomTSS_E14_seq.l)) {
    hm_broad.l[[i]] <- PatternHeatmap(randomTSS_E14_seq.l[[i]], pattern = pattern, coords = range, label = paste(pattern, ", ", samples[i], sep = ""))
  }
  names(hm_broad.l) <- names(randomTSS_E14_seq.l)

# convert heatmaps to metaplot data - returns df 
  meta_broad.l <- lapply(hm_broad.l, function(x) dataMetaplot(x, binsize = 2))

# create a list of dataframes (with sharp and broad column)
  random_df.all <- list()

  samples <- names(randomTSS_E14_seq.l)

  for (i in 1:length(samples)) {
    random_df.all[[i]] <- data.frame("x_coord" = meta_broad.l[[i]]$x_coord, 
                                     "broad_occurrence" = meta_broad.l[[i]]$occurrence,
                                     "sample" = rep(samples[i], times = length(meta_broad.l[[1]]$occurrence)),
                                     "class" = rep("random", times = length(meta_broad.l[[1]]$occurrence)))
  }
  names(random_df.all) <- samples

# collapse the list into a dataframe
  random_df <- do.call(rbind, random_df.all)
  rownames(random_df) <- 1:nrow(random_df)

  all_df$class = "domTSS"

# combine random and normal TSS data
  all_combined.gg <- rbind(all_df, random_df)

# plot metaplots using ggplot2
  library(ggplot2)
  p <- ggplot(all_combined.gg) +
    geom_line(lty = 1, aes(x = x_coord, y = broad_occurrence, colour = class), size = 0.6) +
    scale_color_manual(values = c("royalblue", "goldenrod2"), name = NULL) +
    theme_light()   +
    theme(text = element_text(size = 20, colour = "black"),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          legend.position = "right",
          legend.title = element_blank()) +
    guides(colour = guide_legend(reverse=F))+
    geom_vline(xintercept = 0, lty = "dashed", colour = "gray48", size = 0.6) +
    scale_y_continuous(name = "Relative frequency", limits = c(0, 0.4), breaks = c(0, 0.2, 0.4)) +
    scale_x_continuous(name = "Distance to dominant TSS (bp)", limits = c(-100, 300), seq(-100, 300, by = 100))
  
  p <- p + facet_wrap(~ sample, ncol = 4)
  
  pdf("Figures/WW_shifting_all_domTSS_random.pdf", height = 12, width = 12)
    print(p)
  dev.off()

# zoomed metaplots for the inset = only broad pattern occurence
  library(ggplot2)
  selection.gg <- all_combined.gg[all_combined.gg$sample %in% c("E16_5F", "PN6"), ]
  selection.gg$sample <- droplevels(selection.gg$sample)
  
  p <- ggplot(selection.gg) +
    geom_line(lty = 1, aes(x = x_coord, y = broad_occurrence, colour = class), size = 0.6, legend = FALSE) +
    scale_color_manual(values = c("royalblue4", "goldenrod2")) +
    theme_light()   +
    theme(text = element_text(size = 20, colour = "black"),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          legend.position = "right",
          legend.title = element_blank()) +
    guides(colour = guide_legend(reverse=F)) +
    scale_y_continuous(name = "Relative frequency") +
    scale_x_continuous(name = "Distance to dominant TSS (bp)") 
  
  p <- p + facet_wrap(~ sample, ncol = 2)

  pdf("Figures/WW_shifting_random_sel.pdf", height = 4, width = 11)
    print(p)
  dev.off()

# zoomed metaplots for the inset = only broad pattern occurence for selected samples
  library(ggplot2)
  p <- ggplot(selection.gg) +
    geom_line(lty = 1, aes(x = x_coord, y = broad_occurrence, colour = class), size = 0.6, legend = FALSE) +
    scale_color_manual(values = c("royalblue4", "goldenrod2"), name = NULL) +
    theme_light() +
    theme(text = element_text(size = 20, colour = "black"),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          legend.position = "right",
          legend.title = element_blank()) +
    guides(colour = guide_legend(reverse=F)) +
    scale_y_continuous(name = "Relative frequency") +
    scale_x_continuous(name = "Distance to dominant TSS (bp)") +
    coord_cartesian(xlim = c(50, 200), ylim = c(0.075, 0.15))
  p <- p + facet_wrap(~ sample, ncol = 2)
  
  pdf("Figures/WW_shifting_random_sel.pdf", height = 4, width = 11)
    print(p)
  dev.off()
