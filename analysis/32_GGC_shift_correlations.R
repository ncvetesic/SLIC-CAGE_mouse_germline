#/*==========================================================================#*/
#' ## Extended Data Figure 13E
#/*==========================================================================#*/
# Correlation of GGC distance shifts


## Calculate q0.1 and q0.9 shifts GGCs vc mESC
# import shifting ranges
  PGC_oocyte_embryo_vs_mESC_shift.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/PGC_oocyte_embryo_vs_mESC_shift_grl.RDS")

# import sample specific tag clusters
  tc_PGCs_oocyte_embryo_anno.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/tc_PGCs_oocyte_embryo_anno_grl.RDS")

# extract PGC E16_5F and mESC
  E16_5F_shifts_vs_mESC_shifts.gr <- PGC_oocyte_embryo_vs_mESC_shift.grl[["E16_5F"]]
  E16_5F_TCs.gr <- tc_PGCs_oocyte_embryo_anno.grl[["E16_5F"]]
  mESC_TCs.gr <- tc_PGCs_oocyte_embryo_anno.grl[["E14_mESC"]]

# overlap shifting with sample specific tag cluster - E16_5F
  overlap_E16_5F.pairs <- findOverlapPairs(E16_5F_shifts_vs_mESC_shifts.gr, E16_5F_TCs.gr)
  E16_5F_shifts.gr <- overlap_E16_5F.pairs@first
  mcols(E16_5F_shifts.gr) <- cbind(mcols(E16_5F_shifts.gr), 
                                   "E16_5F_q0.1" =  mcols(overlap_E16_5F.pairs@second)$q_0.1, 
                                   "E16_5F_q0.9" =  mcols(overlap_E16_5F.pairs@second)$q_0.9, 
                                   "E16_5F_tpm" = mcols(overlap_E16_5F.pairs@second)$tpm, 
                                   "E16_5F_IQwidth" = mcols(overlap_E16_5F.pairs@second)$interquantile_width)

# overlap shifting with sample specific tag cluster - mESC
  overlap_mESC.pairs <- findOverlapPairs(E16_5F_shifts.gr, mESC_TCs.gr)
  E16_5F_mESC_shifts.gr <- overlap_mESC.pairs@first
  mcols(E16_5F_mESC_shifts.gr) <- cbind(mcols(E16_5F_mESC_shifts.gr), 
                                        "mESC_q0.1" =  mcols(overlap_mESC.pairs@second)$q_0.1, 
                                        "mESC_q0.9" =  mcols(overlap_mESC.pairs@second)$q_0.9, 
                                        "mESC_tpm" = mcols(overlap_mESC.pairs@second)$tpm, 
                                        "mESC_IQwidth" = mcols(overlap_mESC.pairs@second)$interquantile_width)

# select highest tpm tag clusters in E16.5F - as I have duplicated object
  library(dplyr)
  df <- as.data.frame(E16_5F_mESC_shifts.gr) %>% group_by(start) %>% dplyr::slice(which.max(E16_5F_tpm))
  E16_5F_mESC_shifts_unique.gr <- makeGRangesFromDataFrame(df, keep.extra.columns = T)
  seqinfo(E16_5F_mESC_shifts_unique.gr) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(E16_5F_mESC_shifts_unique.gr)]

# drop chrM
  E16_5F_mESC_shifts_unique.gr <- dropSeqlevels(value = "chrM", E16_5F_mESC_shifts_unique.gr, pruning.mode = "coarse")

# calculate q0.1 distances
  E16_5F_mESC_shifts_unique.gr$q01_dist <- E16_5F_mESC_shifts_unique.gr$E16_5F_q0.1 - E16_5F_mESC_shifts_unique.gr$mESC_q0.1
  E16_5F_mESC_shifts_unique.gr$q01_dist[as.character(strand(E16_5F_mESC_shifts_unique.gr)) == "-"] <- -E16_5F_mESC_shifts_unique.gr$q01_dist[as.character(strand(E16_5F_mESC_shifts_unique.gr)) == "-"]


# calculate q0.9 distances
  E16_5F_mESC_shifts_unique.gr$q09_dist <- E16_5F_mESC_shifts_unique.gr$E16_5F_q0.9 - E16_5F_mESC_shifts_unique.gr$mESC_q0.9
  E16_5F_mESC_shifts_unique.gr$q09_dist[as.character(strand(E16_5F_mESC_shifts_unique.gr)) == "-"] <- -E16_5F_mESC_shifts_unique.gr$q09_dist[as.character(strand(E16_5F_mESC_shifts_unique.gr)) == "-"]

## check if the calculations are correct manually for plus an minus strand - and plus strand, but this kinda shows that they rarely spread in both direcionts either move or spread - goes in favour of 2 alternative nucleosome positions
# later check with double shift and maternal, is the maternal shift usually in the same spot
  plot(E16_5F_mESC_shifts_unique.gr$q01_dist, E16_5F_mESC_shifts_unique.gr$q09_dist, xlim = c(-100, 100), ylim = c(-100, 100))


# have to flip q0.1 and q0.9 on negative strands as for granges q0.1 has to be a lower number, but here I am checking directions so I want the q0.1 to be the upstream range
# extract for pos strand distances
  pos.strands.df <- data.frame("q0.1_dist" = E16_5F_mESC_shifts_unique.gr[as.character(strand(E16_5F_mESC_shifts_unique.gr)) == "+"]$q01_dist,
                               "q0.9_dist" = E16_5F_mESC_shifts_unique.gr[as.character(strand(E16_5F_mESC_shifts_unique.gr)) == "+"]$q09_dist)
  min.strands.df <- data.frame("q0.1_dist" = E16_5F_mESC_shifts_unique.gr[as.character(strand(E16_5F_mESC_shifts_unique.gr)) == "-"]$q09_dist,
                               "q0.9_dist" = E16_5F_mESC_shifts_unique.gr[as.character(strand(E16_5F_mESC_shifts_unique.gr)) == "-"]$q01_dist)

  dist_q01_q09.df <- rbind(pos.strands.df, min.strands.df)

  rm(pos.strands.df)
  rm(min.strands.df)

# plot using ggpplot2
  library(ggplot2)

  p <- ggplot(dist_q01_q09.df, aes(q0.1_dist, q0.9_dist), alpha = 0.3) +
    geom_point(size=3, pch = 21, colour = "black", fill = "black", alpha = 0.3) +
    xlab("q0.1 distance (bp)") +
    ylab("q0.9 distance (bp)") + 
    coord_fixed() +
    theme_bw() +
    theme(text = element_text(size = 16, colour = "black"), 
          axis.text.x = element_text(size = 16, colour = "black"),
          axis.text.y = element_text(size = 16, colour = "black"),
          legend.title = element_blank(),
          legend.position = "none") +
    coord_cartesian(xlim = c(-750, 1250), ylim = c(-750, 1250))


# print plot
  pdf("Figures/quantile_dist_all.pdf", height = 4, width = 4)
    print(p)
  dev.off()

# zoom into the plot
  p <- ggplot(dist_q01_q09.df, aes(q0.1_dist, q0.9_dist), alpha = 0.3) +
    geom_point(size=3, pch = 21, colour = "black", fill = "black", alpha = 0.3) +
    xlab("q0.1 distance (bp)") +
    ylab("q0.9 distance (bp)") + 
    coord_fixed() +
    theme_bw() +
    theme(text = element_text(size = 16, colour = "black"), 
          axis.text.x = element_text(size = 16, colour = "black"),
          axis.text.y = element_text(size = 16, colour = "black"),
          legend.title = element_blank(),
          legend.position = "none") +
    coord_cartesian(xlim = c(-200, 200), ylim = c(-200, 200))
  

# print plot
  pdf("Figures/quantile_dist_all_zoom.pdf", height = 4, width = 4)
    print(p)
  dev.off()


# --------- add spread annotation - upstream_both, downstream_both, up&down, upstream_0, 0_dowstream --------- #
#upstream both - both q0.1 and q0.9 in E16_5F compared to mESC are negative
#downstream both - both q0.1 and q0.9 in E16_5F compared to mESC are positive
#up&down - q0.1 is negative and q0.9 is positive - so it spread in both directions
#upstream_0 - q0.1 is negative so it spread upstream, while q0.9 is unchanged == 0
#0_downstream - q0.1 is unchanged == 0, while q0.9 is positive
#0_0 - no change in iq-position
#down&up - q0.1 is positive and q0.9 is negative - so it narrowed
#down_0 - q0.1 is positive and q0.9 is unchanged - so it narrowed from 5'side
#0_down - q0.1 unchanged and q0.9 is negative - so it narrowed from 3'side


  E16_5F_mESC_shifts_unique.gr$q01_real_dist <- 0  
  E16_5F_mESC_shifts_unique.gr$q09_real_dist <- 0  
  E16_5F_mESC_shifts_unique.gr[as.character(strand(E16_5F_mESC_shifts_unique.gr)) == "+"]$q01_real_dist <- E16_5F_mESC_shifts_unique.gr[as.character(strand(E16_5F_mESC_shifts_unique.gr)) == "+"]$q01_dist
  E16_5F_mESC_shifts_unique.gr[as.character(strand(E16_5F_mESC_shifts_unique.gr)) == "+"]$q09_real_dist <- E16_5F_mESC_shifts_unique.gr[as.character(strand(E16_5F_mESC_shifts_unique.gr)) == "+"]$q09_dist

# reverse q0.1 and q0.9 distances for negative strands
  E16_5F_mESC_shifts_unique.gr[as.character(strand(E16_5F_mESC_shifts_unique.gr)) == "-"]$q09_real_dist <- E16_5F_mESC_shifts_unique.gr[as.character(strand(E16_5F_mESC_shifts_unique.gr)) == "-"]$q01_dist
  E16_5F_mESC_shifts_unique.gr[as.character(strand(E16_5F_mESC_shifts_unique.gr)) == "-"]$q01_real_dist <- E16_5F_mESC_shifts_unique.gr[as.character(strand(E16_5F_mESC_shifts_unique.gr)) == "-"]$q09_dist

# add annotation
  # add upstream both - both q0.1 and q0.9 in E16_5F compared to mESC are negative
    E16_5F_mESC_shifts_unique.gr$shift_anno[E16_5F_mESC_shifts_unique.gr$q01_real_dist < 0 & E16_5F_mESC_shifts_unique.gr$q09_real_dist < 0] <- "upstream_both"

  # add downstream both - both q0.1 and q0.9 in E16_5F compared to mESC are positive
    E16_5F_mESC_shifts_unique.gr$shift_anno[E16_5F_mESC_shifts_unique.gr$q01_real_dist > 0 & E16_5F_mESC_shifts_unique.gr$q09_real_dist > 0] <- "downstream_both"

  # add up&down - q0.1 is negative and a0.9 is positive - so it spread in both directions
    E16_5F_mESC_shifts_unique.gr$shift_anno[E16_5F_mESC_shifts_unique.gr$q01_real_dist < 0 & E16_5F_mESC_shifts_unique.gr$q09_real_dist > 0] <- "up&down"

  # add upstream_0 - q0.1 is negative so it spread upstream, while q0.9 is unchanged == 0
    E16_5F_mESC_shifts_unique.gr$shift_anno[E16_5F_mESC_shifts_unique.gr$q01_real_dist < 0 & E16_5F_mESC_shifts_unique.gr$q09_real_dist == 0] <- "upstream_0"

  # add 0_downstream - q0.1 is unchanged == 0, while q0.9 is positive
    E16_5F_mESC_shifts_unique.gr$shift_anno[E16_5F_mESC_shifts_unique.gr$q01_real_dist == 0 & E16_5F_mESC_shifts_unique.gr$q09_real_dist > 0] <- "0_downstream"

  # add 0_0 - no change in iq-position
    E16_5F_mESC_shifts_unique.gr$shift_anno[E16_5F_mESC_shifts_unique.gr$q01_real_dist == 0 & E16_5F_mESC_shifts_unique.gr$q09_real_dist == 0] <- "0_0"

  # add down&up - q0.1 is positive and q0.9 is negative - so it narrowed
    E16_5F_mESC_shifts_unique.gr$shift_anno[E16_5F_mESC_shifts_unique.gr$q01_real_dist > 0 & E16_5F_mESC_shifts_unique.gr$q09_real_dist < 0] <- "down&up"

  # add down_0 - q0.1 is positive and q0.9 is unchanged - so it narrowed from 5'side
    E16_5F_mESC_shifts_unique.gr$shift_anno[E16_5F_mESC_shifts_unique.gr$q01_real_dist > 0 & E16_5F_mESC_shifts_unique.gr$q09_real_dist == 0] <- "downstream_0"

  # add 0_down - q0.1 unchanged and q0.9 is positive - so it narrowed from 3'side
    E16_5F_mESC_shifts_unique.gr$shift_anno[E16_5F_mESC_shifts_unique.gr$q01_real_dist == 0 & E16_5F_mESC_shifts_unique.gr$q09_real_dist < 0] <- "0_upstream"

  # save the distance annotated E16.5F granges
    saveRDS(E16_5F_mESC_shifts_unique.gr, "intermediate_data/E16_5F_mESC_shifts_unique_distance_anno.gr")

  # organise distances and annotation into a dataframe  
    dist_q01_q09_anno.df <- data.frame("q01_real_dist" =  E16_5F_mESC_shifts_unique.gr$q01_real_dist, 
                                       "q09_real_dist" =  E16_5F_mESC_shifts_unique.gr$q09_real_dist,
                                       "dist_anno" = E16_5F_mESC_shifts_unique.gr$shift_anno)

    names <- c("upstream_both","up&down", "upstream_0", "0_upstream",  "0_0" , "0_downstream", "downstream_0", "down&up", "downstream_both")
    dist_q01_q09_anno.df$dist_anno <- factor(dist_q01_q09_anno.df$dist_anno, levels = names)

# colorscheme for all options 9 - red to blue diverging so organise from most most upstream to downstream
  col <- c("#b2182b", "#8073ac","#fddbc7","#878787","white", "#abd9e9", "#74add1", "black", "#2166ac")
  names(col) <- names


# zoom into the plot
p <- ggplot(dist_q01_q09_anno.df, aes(q01_real_dist, q09_real_dist, fill = dist_anno), alpha = 0.6) +
  geom_point(size=3, pch = 21, colour = "black", alpha = 0.6) +
  scale_fill_manual(values = col) +
  xlab("q0.1 distance (bp)") +
  ylab("q0.9 distance (bp)") + 
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size = 16, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        legend.title = element_blank(),
        legend.position = "right") +
  coord_cartesian(xlim = c(-200, 200), ylim = c(-200, 200)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", size = 0.7)


# print plot
  pdf("Figures/quantile_dist_all_zoom_anno.pdf", height = 6, width = 9)
    print(p)
  dev.off()

# calculate correlation
  cor(dist_q01_q09_anno.df$q01_real_dist, dist_q01_q09_anno.df$q09_real_dist, method = "pearson")

# pearson corr 0.9345018

  table(dist_q01_q09_anno.df$dist_anno)
