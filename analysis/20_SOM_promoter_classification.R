#/*==========================================================================#*/
#' ## Extended Data Figure 8A
#/*==========================================================================#*/
# SOM promoter classification

## Promoters only - custom SOM plotting on promoters only tag clusters - need consensus clusters annotated
  library(kohonen)
  library(tibble)
  library(stringr)
  library(forcats)
  
# Annotate consensus clusters to select promoter only consensus clusters
  # load  consensus clusters
    consensus.clusters.gr <- readRDS("intermediate_data/consensus.clusters.gr.RDS")
  
  # add clustId information
    consensus.clusters.gr$cons_id <- as.numeric(1:length(consensus.clusters.gr))
  
    library(ChIPseeker)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  
  # use peak anno from ChipSeeker to annotate switching tag clusters
    peakAnno.l <- lapply(GRangesList(consensus.clusters.gr), function(x) annotatePeak(x, TxDb = txdb,  annoDb = NULL, sameStrand = TRUE, verbose = FALSE))
  
  # convert to GRanges 
    consensus.clusters.anno.gr <- lapply(peakAnno.l, as.GRanges)
  
  # save annotated GRanges object 
    saveRDS(consensus.clusters.anno.gr, "intermediate_data/consensus.clusters.anno.gr.RDS")
  
  # plot annotated feature
    png(file = "merged_replicates/results/consensus_clusters_anno.png", height = 250, width = 800)
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

  
  # order by % of each feature
    order_lev <- unique(feats.df$Feature[order(feats.df$Frequency, decreasing = FALSE)])
    feats.df$Feature <- factor(feats.df$Feature, levels = order_lev)
  

  # plot genomic features
    library(ggplot2)
    p <- ggplot(feats.df, aes(x = Feature, y = Frequency, fill = Feature), alpha = 0.7) +
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
    
    pdf(file = "merged_replicates/results/cons_clust_ggplot_gen_loc.pdf", height = 4, width = 6)
      print(p)
    dev.off()

# SOM classification of consensus clusters  
  # load annotated consensus clusters
    consensus.clusters.anno.gr <- readRDS("intermediate_data/consensus.clusters.anno.gr.RDS")
    consensus.clusters.anno.gr <- consensus.clusters.anno.gr[[1]]

  # select promoter regions only
    promoters <- c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")
    consensus.clusters.prom.gr <- consensus.clusters.anno.gr[consensus.clusters.anno.gr$annotation %in% promoters]

# consensus cluster level
  CAGE.sonc.tpm <- consensusClustersTpm(CAGEset_PGC_embryo_merged)

# select promoters only from matrix
  CAGE.conc.prom.tpm <- CAGE.sonc.tpm[consensus.clusters.prom.gr$cons_id, ]

## CAGEr normalization
  norm.for.som <- function(tpm.mx, tpmThreshold = 1, nrPassThreshold = 1){
    sample.labels <- colnames(tpm.mx)
    nr.pass.threshold <- apply(tpm.mx, 1, function(x) {sum(x >= tpmThreshold)})
    idx <- nr.pass.threshold >= min(nrPassThreshold, length(sample.labels))
    ## scale divides the column by the scaling factor
    ## the scaling factor is calculate by applying sqrt(sum(x$ESCday0^2)/(n-1)), where n = number of rows with non-missing values
    ## non-negative log(tpm) are scaled by 
    ## by transposing the matrix with t, the scaling factor is calculated for each element over all the samples
    matForSom <- t(base::scale(t(log(tpm.mx+1)), center=F)) 
    matForSom <- matForSom[idx,]
    return(matForSom)
  }

  motForSom_Cons <- norm.for.som(CAGE.conc.prom.tpm)
  motForSom_Cons <- na.omit(motForSom_Cons)

## som here is needed to cluster each element based on normalised score
  openSom_cons <- kohonen::som(motForSom_Cons,
                               grid = somgrid(3, 3, 'hexagonal'),
                               mode = 'pbatch',
                               cores = 4)

  somDf_cons <- as.data.frame(motForSom_Cons)
  somDf_cons <- rownames_to_column(somDf_cons, var = "elID")
  somDf_cons$som <- openSom_cons$unit.classif

# windsorize signal 
  somDf_cons_wins <- somDf_cons %>% dplyr::group_by(som) %>% dplyr::select(-c(som, elID)) %>% purrr::map(function(x) Winsorize(x)) %>% as.data.frame() 
  somDf_cons_wins$elID <- somDf_cons$elID

## for consensusClusteres
  somDf_cons <- gather(somDf_cons, "sample", "value", -c("elID", "som")) # transposes the columns in rows
  somDf_cons_wins <- gather(somDf_cons_wins, "sample", "value", -c("elID", "som")) # transposes the columns in rows


  somStats_cons <- somDf_cons %>% group_by(som, sample) %>% dplyr::summarize(avg = median(value))
  somStats_cons_wins <- somDf_cons_wins %>% group_by(som, sample) %>% dplyr::summarize(avg = median(value))
  
  somStats_cons$sample <- factor(somStats_cons$sample, levels = sampleLabels(CAGEset_PGC_embryo_merged))
  somStats_cons_wins$sample <- factor(somStats_cons_wins$sample, levels = sampleLabels(CAGEset_PGC_embryo_merged))


  somDf_cons$sample <- factor(somDf_cons$sample, levels = sampleLabels(CAGEset_PGC_embryo_merged))
  somDf_cons_wins$sample <- factor(somDf_cons_wins$sample, levels = sampleLabels(CAGEset_PGC_embryo_merged))

# count the number of clusters per SOM
  som_count <- somDf_cons %>% dplyr::group_by(som, sample) %>% dplyr::select(elID) %>% dplyr::tally()
  som_count_summary <- data.frame("som "= unique(som_count$som), "number" = unique(som_count$n))

  som_count_wins <- somDf_cons_wins %>% dplyr::group_by(som, sample) %>% dplyr::select(elID) %>% dplyr::tally()
  som_count_summary_wins <- data.frame("som "= unique(som_count_wins$som), "number" = unique(som_count_wins$n))

# add number to som class
  somDf_cons$som <- paste0(somDf_cons$som, " (", som_count_summary[somDf_cons$som, ]$number, ")")
  somStats_cons$som <- paste0(somStats_cons$som, " (", som_count_summary[somStats_cons$som, ]$number, ")")

  somDf_cons_wins$som <- paste0(somDf_cons_wins$som, " (", som_count_summary_wins[somDf_cons_wins$som, ]$number, ")")
  somStats_cons_wins$som <- paste0(somStats_cons_wins$som, " (", som_count_summary_wins[somStats_cons_wins$som, ]$number, ")")

##plot
  p <-  ggplot(somDf_cons_wins, aes(x = sample, y = value, fill = sample), inherit.aes = FALSE) + 
    geom_violin(alpha = .8, trim = TRUE, lwd = 0.25) +
    scale_fill_manual(values = c("firebrick1",  "#FDE725FF",  "#E8E419FF",  "#D1E11CFF",  "#B9DE28FF",  "#A2DA37FF",  "#8BD646FF",  "#75D054FF",
                                 "#61CB5FFF",  "#4FC46AFF",  "#3FBC73FF",  "#31B57BFF",  "#440154FF", "#46337EFF",  "#3F4889FF", 
                                 "#365C8DFF",  "#2E6E8EFF" )) +
    facet_wrap(~ som) + 
    geom_point(data = somStats_cons_wins, aes(x = sample, y = avg), inherit.aes = FALSE, size = 0.5) + 
    geom_line(data = somStats_cons_wins, aes(group = 1, x = sample, y = avg), inherit.aes = FALSE, lwd = 0.25) +
    scale_color_manual(values = col_scheme) +
    theme_bw() +
    theme(text = element_text(size = 16),
          axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 1, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          strip.background =element_rect(fill="white")) +
    coord_cartesian(ylim=(c(0, 5)))
  
  pdf("Figures/SOM_samples_promOnly.pdf", height = 8, width = 12)
  print(p)
  dev.off() 


# ----- extract SOM genes ----- #
  # consensus cluster ID corresponds to elID
  ## extract consensus clusters
    consensus.clusters <- consensusClusters(CAGEset_PGC_embryo_merged)
  
  # convert consensus clusters into GRanges
    consensus.clusters.gr <- GRanges(seqnames = consensus.clusters$chr,
                                     range = IRanges(start = consensus.clusters$start, end = consensus.clusters$end),
                                     strand = consensus.clusters$strand,
                                     tpm = consensus.clusters$tpm)
  
  # add genome info
    seqinfo(consensus.clusters.gr) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(consensus.clusters.gr)]
  
  # split consensus cluster IDs according to SOM class
  # soms
    som_df.l <- split(somDf_cons_wins, list(somDf_cons_wins$som, somDf_cons_wins$sample))
  
  # extract consensus cluster IDs - all are the same regardless of sample, so I will extract 9 SOM ids from mESC
    consID_som.l <- list()
    for(i in 1:9) {
      consID_som.l[[i]] <- as.numeric(som_df.l[[i]]$elID)
    }

# divide consensus clusters according to ID
  consensus.clusters.som.l <- list()

  for(i in 1:length(consID_som.l)) {
    consensus.clusters.som.l[[i]] <- consensus.clusters.gr[consID_som.l[[i]], ]
  }

# convert go GRangesList
  consensus.clusters.som_promOnly.grl <- GRangesList(consensus.clusters.som.l)
  names(consensus.clusters.som_promOnly.grl) <- paste0("som", 1:9)

# save SOM classes - Damir code
  saveRDS(consensus.clusters.som_promOnly.grl, "intermediate_data/consensus_clusters_som_promOnly_Damir.9classes_grl.RDS")
```