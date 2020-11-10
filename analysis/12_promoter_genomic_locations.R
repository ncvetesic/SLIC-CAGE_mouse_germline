#/*==========================================================================#*/
#' ## Extended Data Figure 2D
#/*==========================================================================#*/
# genomic locations of tag clusters identified in each sample (tag clusters, i.e. promoters)

# load libraries
  library(ChIPseeker)
  library(GenomicFeatures)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# create handle for txdb mouse
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# load tag clusters
  tc_PGCs_oocyte_embryo.grl <- readRDS(file ="../all_reps_PN6_2cell/intermediate_data/tc_PGCs_oocyte_embryo_grl.RDS")

# use peak anno from ChipSeeker to annotate tag clusters
  peakAnno.l <- lapply(tc_PGCs_oocyte_embryo.grl, function(x) {
    annotatePeak(x, TxDb = txdb,  annoDb = NULL, tssRegion = c(-500, 100), sameStrand = TRUE, verbose = FALSE)
  })

# convert to granges
  tc_PGCs_oocyte_embryo_anno.grl <- lapply(peakAnno.l, as.GRanges)

# add sequence information
  for(i in 1:length(tc_PGCs_oocyte_embryo_anno.grl)) {
    seqinfo(tc_PGCs_oocyte_embryo_anno.grl[[i]]) <- seqinfo(tc_PGCs_oocyte_embryo_anno.grl[[i]])
  }

# plot annotated feature
  png(file = "../all_reps_PN6_2cell/genomic_location/annoFeat.png", height = 500, width = 800)
  plotAnnoBar(peakAnno.l, title = NULL)
  dev.off()

# plotting is limited within the package so I extract genomic locations
  feats.l <- lapply(peakAnno.l, function(x) {
  return(x@annoStat)
})

# round percentages
  feats.l <- lapply(feats.l, function(x) {
    x$Frequency <- round(as.numeric(x$Frequency), digits = 2)
    return(x) 
  })

# new colour scheme
  col <- c("#F8F2AB", "#B4D6A4", "#67BC9A", "#13B0A5", "#0071A7", "#3E606F", "#88F9D4", "#18C29C",
           "#0B877D", "#126872", "#031727")

# set feature factors
  feature_factors <- c("Promoter", "Promoter (<=1kb)", "Promoter (1-3kb)", "5' UTR", "1st Exon", "Other Exon", 
                       "1st Intron", "Other Intron", "3' UTR", "Downstream (<=300)", "Distal Intergenic")
  names(col) <- feature_factors
  col_sel <- col[names(col) %in% unique(feats.l[[1]]$Feature)]

# order by % of each feature
  for(i in 1:length(feats.l)) {
    order_lev <- feats.l[[i]]$Feature[order(feats.l[[i]]$Frequency, decreasing = FALSE)]
    feats.l[[i]]$Feature <- factor(feats.l[[i]]$Feature, levels = order_lev)
  }

# plot genomic features
  library(ggplot2)
  for(i in 1:length(feats.l)) {
    p <- ggplot(feats.l[[i]], aes(x = Feature, y = Frequency, fill = Feature), alpha = 0.7) +
      geom_bar(stat = "identity", width = 0.75, colour = "black", lwd = 0.125) +
      coord_flip() +
      scale_fill_manual("Features", values = col_sel) +
      theme_bw() +
      theme(text = element_text(size = 14, colour = "black"),
            legend.title = element_blank(),
            axis.title.x = element_text(colour = "black"),
            axis.title.y = element_text(colour = "black"),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      labs(y = "Percentage", x = NULL) +
      guides(fill = guide_legend(reverse = TRUE)) +
      ylim(0, 100)
  
  pdf(file = paste0("../all_reps_PN6_2cell/genomic_location/tc_genomic_loc_reanno_CAGE_", names(feats.l)[i], ".pdf"), height = 3, width = 5)
  print(p)
  dev.off()
  }

# plot genomic features all together
# add sample name
  feats.l <- lapply(1:length(feats.l), function(x) {
    df <- feats.l[[x]]
    df$sample <- names(feats.l)[x]
    return(df)
  })
  names(feats.l) <- names(tc_PGCs_oocyte_embryo_anno.grl)

# flatten list to a single dataframe
  feats.df <- do.call(rbind, feats.l)
  feats.df$sample <- factor(feats.df$sample, levels = rev(names(tc_PGCs_oocyte_embryo_anno.grl)))

# print genomic features
  library(ggplot2)
  p <- ggplot(feats.df, aes(x = sample, y = Frequency, fill = Feature), alpha = 0.7) +
    geom_bar(stat = "identity", width = 0.75, colour = "black", lwd = 0.125) +
    coord_flip() +
    scale_fill_manual("Features", values = col_sel) +
    theme_bw() +
    theme(text = element_text(size = 14, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(y = "Percentage", x = NULL) +
    guides(fill = guide_legend(reverse = TRUE)) +
    ylim(0, 101)

pdf(file = "../April_2019_replicates_oocyte/Figures/tc_genomic_loc_mESC_PGC_oocyte.pdf", height = 4, width = 6)
  print(p)
dev.off()