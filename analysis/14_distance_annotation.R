#/*==========================================================================#*/
#' ## Extended Data Figure 3A-D
#/*==========================================================================#*/
# distance of identified oocyte tag cluster from mm10 annotation or oocyte annotation from Veselovska et al 2015 (Genome Biol 16, 209)

# its grcm38 which is mm10
# should reannotate using chippeakanno - upstream of mm10 or upstream of GK gtfs
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(rtracklayer)
  oocyte_Kelsey_gtf.gr <- import("../../extra_data/13059_2015_769_MOESM5_ESM.gtf") #available in the publication - see above

  seqlevels(oocyte_Kelsey_gtf.gr) <- paste0("chr", seqlevels(oocyte_Kelsey_gtf.gr))
  seqinfo(oocyte_Kelsey_gtf.gr) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(oocyte_Kelsey_gtf.gr)]

# select transcripts
  oocyte_Kelsey_gtf_transcripts.gr <- oocyte_Kelsey_gtf.gr[oocyte_Kelsey_gtf.gr$type %in% "transcript"]
  transcripts_mm10.gr <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)

# annotate oocyte samples with Kelsey annotation
# read tag clusters
  tc_PGCs_oocyte_embryo.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/tc_PGCs_oocyte_embryo_grl.RDS")

  tc_oocytes.grl <- tc_PGCs_oocyte_embryo.grl[c("PN6", "PN14", "oocyte")]

# annotate using chippeakanno - note the txdb is gavin kelsey gtf based - and TSSs should be overlapping or upstream
  library(ChIPpeakAnno)
  tc_oocytes_GK_anno.grl <- lapply(tc_oocytes.grl, function(x) {
    annotatePeakInBatch(x, AnnotationData = oocyte_Kelsey_gtf_transcripts.gr, output = c("nearestLocation"),
                        ignore.strand = F, FeatureLocForDistance="TSS", bindingRegion = c(-500, 100), multiple = F)
  })

# change annotation colnames to distinguish
  tc_oocytes_GK_anno_rc.grl <- lapply(tc_oocytes_GK_anno.grl, function(x){
    colnames(mcols(x)) <- c("nr_ctss", "dominant_ctss", "tpm",  "tpm.dominant_ctss", "q_0.1", "q_0.9",  "interquantile_width",  "GK_peak", "GK_feature", "GK_start_position", "GK_end_position", "GK_feature_strand", "GK_insideFeature", "GK_distancetoFeature", "GK_shortestDistance", "GK_fromOverlappingOrNearest")
    return(x)
  })


# annotate using chippeakanno - note the txdb is gavin kelsey gtf based - and TSSs should be overlapping or upstream
  library(ChIPpeakAnno)
  tc_oocytes_GK_mm10_anno.grl <- lapply(tc_oocytes_GK_anno_rc.grl, function(x) {
    annotatePeakInBatch(x, AnnotationData = transcripts_mm10.gr, output = c("nearestLocation"),
                        ignore.strand = F, FeatureLocForDistance="TSS", bindingRegion = c(-500, 100), multiple = F)
  })


# save annotated GRanges object 
  saveRDS(tc_oocytes_GK_mm10_anno.grl, "merged_replicates/intermediate_data/tc_oocytes_GK_mm10_anno_grl.RDS")

# plot annotated feature
# --- GK 
  pdf(file = "Figures/pie_GK_annotation_p6.pdf", height = 5, width = 5)
  pie1(table(as.data.frame(tc_oocytes_GK_mm10_anno.grl[[1]]$GK_insideFeature))) 
  dev.off()

  pdf(file = "Figures/pie_GK_annotation_p14.pdf", height = 5, width = 5)
  pie1(table(as.data.frame(tc_oocytes_GK_mm10_anno.grl[[2]]$GK_insideFeature))) 
  dev.off()

  pdf(file = "Figures/pie_GK_annotation_Mii.pdf", height = 5, width = 5)
  pie1(table(as.data.frame(tc_oocytes_GK_mm10_anno.grl[[3]]$GK_insideFeature))) 
  dev.off()
  
# --- mm10 
  pdf(file = "Figures/pie_mm10_annotation_p6.pdf", height = 5, width = 5)
  pie1(table(as.data.frame(tc_oocytes_GK_mm10_anno.grl[[1]]$insideFeature))) 
  dev.off()
  
  pdf(file = "Figures/pie_mm10_annotation_p14.pdf", height = 5, width = 5)
  pie1(table(as.data.frame(tc_oocytes_GK_mm10_anno.grl[[2]]$insideFeature))) 
  dev.off()
  
  pdf(file = "Figures/pie_mm10_annotation_Mii.pdf", height = 5, width = 5)
  pie1(table(as.data.frame(tc_oocytes_GK_mm10_anno.grl[[3]]$insideFeature))) 
  dev.off()

# extract distances
  distances.l <- lapply(1:length(tc_oocytes_GK_mm10_anno.grl), function(x) {
    distance.df <- data.frame(distanceToFeature_GK = tc_oocytes_GK_mm10_anno.grl[[x]]$GK_distancetoFeature,
                              distanceToFeature_mm10 = tc_oocytes_GK_mm10_anno.grl[[x]]$distancetoFeature,
                              sample = names(tc_oocytes_GK_mm10_anno.grl)[x])
    return(distance.df)
  })

# flatten lists to dataframes
  distances.df <- do.call(rbind, distances.l)

# change sample level
  distances.df$sample <- factor(distances.df$sample, levels = c("PN6", "PN14", "oocyte"))

# flatten
  library(tidyr)
  distances.gg <- gather(distances.df, "type", "distance", -c("sample"))
  distances.gg$type <- factor(distances.gg$type, levels = c("distanceToFeature_GK", "distanceToFeature_mm10"))

# remove NA
  distances_filtr.gg <- distances.gg[!is.na(distances.gg$distance), ]

# convert to abs values
  distances_filtr.gg$distance <- abs(distances_filtr.gg$distance)

# plot  boxplot of distances
  library(ggplot2)
  p <- ggplot(distances_filtr.gg, aes(y = log2(distance +1), x = sample , fill = type)) +
    geom_boxplot(outlier.shape = NA, varwidth = T, notch = F, alpha = 0.6) +
    scale_fill_manual(values = c("gray66", "cornflowerblue")) +
    scale_x_discrete(labels = levels(distances_filtr.gg$sample)) +
    theme_bw() +
    theme(text = element_text(size = 16),
          axis.text.x = element_text(size = 16, angle = 45, hjust = 1, colour = "black"),
          axis.text.y = element_text(size = 16,  colour = "black"),
          axis.title.x = element_text(size = 16,  colour = "black"),
          axis.title.y = element_text(size = 16,  colour = "black"),
          legend.position = "right") +
    labs(title = "") +
    labs( x = "", y = "log2(abs(distance) + 1)") +
    geom_hline(yintercept = median(log2(distances_filtr.gg$distance + 1)), linetype="dashed", color = "indianred1", size=1, alpha = 0.8)
  
  pdf("Figures/boxplot_annot_dist.pdf", height = 5, width = 6)
    print(p)
  dev.off() 

# plot heatmaps coverage from transcripts = create domTSS heatmap position from mm10 and GK transcripts
  oocyte_Kelsey_gtf_transcripts_TSS.gr <- promoters(oocyte_Kelsey_gtf_transcripts.gr, upstream = 0, downstream = 1)
  transcripts_mm10_TSS.gr <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)

# create centring from GK annotated object so I can use distance for sorting
# filter NA out
  tc_oocytes_GK_mm10_anno_filtr.grl <- lapply(tc_oocytes_GK_mm10_anno.grl, function(x) {
    gr <- x[!is.na(x$GK_feature)]
    return(gr)
  })

  tc_oocytes_GK_TSS.grl <- lapply(tc_oocytes_GK_mm10_anno_filtr.grl, function(x){
    gr <- GRanges(seqnames = seqnames(x),
                  ranges = IRanges(start = mcols(x)$GK_start_position, 
                                   end = mcols(x)$GK_end_position),
                  strand = strand(x))
    mcols(gr) <- mcols(x)
  
  prom.gr <- promoters(gr, upstream = 0, downstream = 1)
  return(prom.gr)
  })

# centring from mm10
  tc_oocytes_mm10_anno_filtr.grl <- lapply(tc_oocytes_GK_mm10_anno.grl, function(x) {
    gr <- x[!is.na(x$feature)]
    return(gr)
  })

  tc_oocytes_mm10_TSS.grl <- lapply(tc_oocytes_mm10_anno_filtr.grl, function(x){
    gr <- GRanges(seqnames = seqnames(x),
                  ranges = IRanges(start = mcols(x)$start_position, 
                                   end = mcols(x)$end_position),
                  strand = strand(x))
    mcols(gr) <- mcols(x)
    
    prom.gr <- promoters(gr, upstream = 0, downstream = 1)
    return(prom.gr)
  })

# add seqinfo and drop chrM
  for(i in 1:length(tc_oocytes_mm10_TSS.grl)) {
    seqinfo(tc_oocytes_mm10_TSS.grl[[i]]) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(tc_oocytes_mm10_TSS.grl[[i]])]
    seqinfo(tc_oocytes_GK_TSS.grl[[i]]) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(tc_oocytes_GK_TSS.grl[[i]])]
    tc_oocytes_mm10_TSS.grl[[i]] <- dropSeqlevels(tc_oocytes_mm10_TSS.grl[[i]], "chrM", pruning.mode = "coarse")
    tc_oocytes_GK_TSS.grl[[i]] <- dropSeqlevels(tc_oocytes_GK_TSS.grl[[i]], "chrM", pruning.mode = "coarse")
    tc_oocytes.grl[[i]] <- dropSeqlevels(tc_oocytes.grl[[i]], "chrM", pruning.mode = "coarse")
  }


## plot heatmaps centred on mm10 start transcripts and order by distance from the starts
# ---- centring from mm10
# create windows for plotting
  up <- 1000
  down <- 1000
  range <- c(-up, down)
  win <- up + down

  domTSS_win.grl <- lapply(tc_oocytes_mm10_TSS.grl, function(x) promoters(x, up = up, down = down))
  
  # remove out of bound ranges
  domTSS_win.grl <- lapply(domTSS_win.grl, function(x) x[width(trim(x)) == win])
  
  # check the length of the each range in the list
  sapply(domTSS_win.grl, length)
  
  # randomize prior to sorting
  set.seed(10)
  random_row.l <- lapply(sapply(domTSS_win.grl, length), sample)
  
  for (i in 1:length(domTSS_win.grl)) {
    domTSS_win.grl[[i]] <- domTSS_win.grl[[i]][random_row.l[[i]]]  
  }
  
  # sort index - no sequences are trimmed out so I can use the .gr object for sorting - using distance to feature
  sort_idx.l <- lapply(domTSS_win.grl, function(x) order(x$distancetoFeature,  decreasing = FALSE))
  for (i in 1:length(domTSS_win.grl)) {
    domTSS_win.grl[[i]] <- domTSS_win.grl[[i]][sort_idx.l[[i]]]  
  }
  
  for (i in 1:length(domTSS_win.grl))  {
    hm.l = CoverageHeatmap(
      domTSS_win.grl[[i]],
      tc_oocytes.grl[[i]],
      coords = range,
      weight = log2(tc_oocytes.grl[[i]]$tpm + 1))
    hm_smoothed.l <- smoothHeatmap(hm.l, sigma = c(1, 1), output.size=c(20000, 1000))
    
    scale(hm_smoothed.l) <- c(0, 5)
    
    # plot heatmap
    pdf(paste0("Figures/GK_annotation/mm10_tc_cov_dist_to_feat_", names(domTSS_win.grl)[i], ".pdf"), height = 5.5, width = 4)
    plotHeatmapList(hm_smoothed.l,
                    legend = F,
                    color = "Greys")
    dev.off()
  }


## plot heatmaps centred on GK start transcripts and order by distance from the starts
# ---- centring from mm10
# create windows for plotting
  up <- 1000
  down <- 1000
  range <- c(-up, down)
  win <- up + down
  
  domTSS_win.grl <- lapply(tc_oocytes_GK_TSS.grl, function(x) promoters(x, up = up, down = down))
  
  # remove out of bound ranges
  domTSS_win.grl <- lapply(domTSS_win.grl, function(x) x[width(trim(x)) == win])
  
  # check the length of the each range in the list
  sapply(domTSS_win.grl, length)
  
  # randomize prior to sorting
  set.seed(10)
  random_row.l <- lapply(sapply(domTSS_win.grl, length), sample)
  
  for (i in 1:length(domTSS_win.grl)) {
    domTSS_win.grl[[i]] <- domTSS_win.grl[[i]][random_row.l[[i]]]  
  }
  
  # sort index - no sequences are trimmed out so I can use the .gr object for sorting - using distance to feature
  sort_idx.l <- lapply(domTSS_win.grl, function(x) order(x$GK_distancetoFeature,  decreasing = FALSE))
  for (i in 1:length(domTSS_win.grl)) {
    domTSS_win.grl[[i]] <- domTSS_win.grl[[i]][sort_idx.l[[i]]]  
  }
  
  for (i in 1:length(domTSS_win.grl))  {
    hm.l = CoverageHeatmap(
      domTSS_win.grl[[i]],
      tc_oocytes.grl[[i]],
      coords = range,
      weight = log2(tc_oocytes.grl[[i]]$tpm + 1))
    hm_smoothed.l <- smoothHeatmap(hm.l, sigma = c(1, 1), output.size=c(20000, 1000))
    scale(hm_smoothed.l) <- c(0, 5)
    
    # plot heatmap
    pdf(paste0("Figures/GK_annotation/GK_tc_cov_dist_to_feat_", names(domTSS_win.grl)[i], ".pdf"), height = 5.5, width = 4)
    plotHeatmapList(hm_smoothed.l,
                    legend = F,
                    color = "Greys")
    dev.off()
  }


# plot feature distance distributions to show asymetry in P6 sample; objects needed: tc_oocytes_GK_mm10_anno_filtr.grl, tc_oocytes_mm10_TSS.grl
# calculate distances again because above i used abs values for boxplot
# flatten lists to dataframes
  distances.df <- do.call(rbind, distances.l)

# change sample level
  distances.df$sample <- factor(distances.df$sample, levels = c("PN6", "PN14", "oocyte"))

# flatten
  library(tidyr)
  distances.gg <- gather(distances.df, "type", "distance", -c("sample"))
  distances.gg$type <- factor(distances.gg$type, levels = c("distanceToFeature_GK", "distanceToFeature_mm10"))

# remove NA
  distances_filtr.gg <- distances.gg[!is.na(distances.gg$distance), ]

library(dplyr)  
  P6_distances.gg <- distances_filtr.gg %>% filter(sample == "PN6")
  P6_distances.gg$sample <- droplevels(P6_distances.gg$sample)
  
  # calculate medians
  P6_dist_GK.med <- P6_distances.gg %>% filter(type == "distanceToFeature_GK") 
  P6_dist_GK.med <- median(P6_dist_GK.med$distance)
  
  P6_dist_mm10.med <- P6_distances.gg %>% filter(type == "distanceToFeature_mm10")
  P6_dist_mm10.med <- median(P6_dist_mm10.med$distance)
  
  # print plots
  library(ggplot2)
  
  p <- ggplot(P6_distances.gg, aes(x = distance,  fill = type), alpha = 0.2) +
    geom_histogram(aes(x = distance, y = ..count../sum(..count..)*100), binwidth = 5, col = "black", lwd = 0.5, position = "identity", alpha = 0.4) +
    scale_fill_manual("type", values = c("gray66", "cornflowerblue")) +
    theme_bw()  +
    theme(text = element_text(size = 16, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          legend.position = "top") +
    labs(y = "%", x = "Distance (bp)") + 
    coord_cartesian( xlim = c(-250, 250), ylim = c(0, 2.5)) +
    geom_vline(xintercept = P6_dist_GK.med, lty = "dashed", lwd = 1.25, col = "gray66") +
    geom_vline(xintercept = P6_dist_mm10.med, lty = "dashed", lwd = 1.25, col = "cornflowerblue")  

  
  pdf("Figures/Kelsey_mm10_anno_P6_dist_hist.pdf", width = 7, height = 4)
    print(p) 
  dev.off() 