#/*==========================================================================#*/
#' ## Extended Data Figure 8B
#/*==========================================================================#*/
# Genomic locations of SOM promoter classes

### Annotation of SOM classes
# load som consensus clusters
  consensus.clusters.som_promOnly.grl <- readRDS("intermediate_data/consensus_clusters_som_promOnly_Damir.9classes_grl.RDS")

  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# use peak anno from ChipSeeker to annotate switching tag clusters
  peakAnno.l <- lapply(GRangesList(consensus.clusters.som_promOnly.grl), function(x) annotatePeak(x, TxDb = txdb,  annoDb = NULL, sameStrand = TRUE, verbose = FALSE))

# convert to GRanges 
  consensus.clusters.som_promOnly.anno.grl <- lapply(peakAnno.l, as.GRanges)

# save annotated GRanges object 
  saveRDS(consensus.clusters.som_promOnly.anno.grl, "intermediate_data/consensus_clusters_som_promOnly_anno_Damir.9classes_grl.RDS")

# plot annotated feature
  png(file = "merged_replicates/results/SOM/consensus_clusters_SOM9_promOnly_Damir.png", height = 800, width = 800)
    plotAnnoBar(peakAnno.l, title = NULL)
  dev.off()

# plotting is limited within the package so I extract genomic locations
  feats.l <- lapply(peakAnno.l, function(x) return(x@annoStat))

# set class
  for (i in 1:length(feats.l)) {
    feats.l[[i]]$Frequency <- round(as.numeric(feats.l[[i]]$Frequency), digits = 2)
  }

# colour scheme
  col <- c("#F8F2AB","#FF9F27", "#B4D6A4", "#67BC9A", "#13B0A5", "#0071A7", "#3E606F", "#88F9D4", "#18C29C",
           "#0B877D", "#126872", "#031727")

# set feature factors
  feature_factors <- c("Promoter", "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR", "1st Exon", "Other Exon", "1st Intron", "Other Intron", "3' UTR", "Downstream (<=3kb)", "Distal Intergenic")

  names(col) <- feature_factors
  col_sel <- col[names(col) %in% unique(unlist(lapply(feats.l, function(x) x$Feature)))]

# convert the list to a dataframe
  feats.df <- do.call(rbind, feats.l)

# add column of sample names
  soms <- sapply(strsplit(rownames(feats.df), "\\." ), function(x) return(x[1]))
  feats.df$som <- soms

# order by % of each feature
  order_lev <- unique(feats.df$Feature[order(feats.df$Frequency, decreasing = FALSE)])
  feats.df$Feature <- factor(feats.df$Feature, levels = order_lev)

# set levels for samples
  feats.df$som <- factor(feats.df$som, levels = paste0("som", 9:1))

# plot genomic features
  library(ggplot2)
  p <- ggplot(feats.df, aes(x = som, y = Frequency, fill = Feature), alpha = 0.7) +
    geom_bar(stat = "identity", width = 0.75, colour = "black", lwd = 0.125) +
    coord_flip() +
    scale_fill_manual("Features", values = col_sel) +
    theme_bw() +
    theme(text = element_text(size = 12, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(y = "Percentage", x = NULL) +
    guides(fill = guide_legend(reverse = TRUE)) +
    ylim(0, 100.1)

  pdf(file = "Figures/cons_clust_SOM9_promOnly_Damir_ggplot_gen_loc.pdf", height = 4, width = 6)
    print(p)
  dev.off()
