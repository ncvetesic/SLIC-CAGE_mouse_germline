#/*==========================================================================#*/
#' ## Figure 2C-F, 3B,C, 4D,F,G
#/*==========================================================================#*/
# Visualisation of sequence patterns using heatmaps - all visualisations are similar, different centring or ordering is used, below are heatmap examples

## PLotting WW, SS and TATAWAWR patterns in shifting promoters
  library(heatmaps)
  library(BSgenome.Mmusculus.UCSC.mm10)
  domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/domTSS_PGC_oocyte_embryo_vs_mESC_shift_grl.RDS")
  domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/domTSS_PGC_oocyte_embryo_vs_mESC_shift_shift_X_grl.RDS")
  PGC_oocyte_embryo_vs_mESC_shift_anno.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/PGC_oocyte_embryo_vs_mESC_shift_anno_grl.RDS")

# keep standard chromosomes  
  domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl <- lapply(domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl, function(x){
    gr <- dropSeqlevels(x, "chrM", pruning.mode = "coarse")
  })
  
  domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl <- lapply(domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl, function(x){
    gr <- dropSeqlevels(x, "chrM", pruning.mode = "coarse")
  })
  
  PGC_oocyte_embryo_vs_mESC_shift_anno.grl <- lapply(PGC_oocyte_embryo_vs_mESC_shift_anno.grl, function(x){
    gr <- dropSeqlevels(x, "chrM", pruning.mode = "coarse")
  })

  samples <- names(domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl)

# create windows for plotting
  up <- 500
  down <- 500
  range <- c(-up, down)
  win <- up + down

# centre on mESC
  domTSS_win.grl <- lapply(domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl, function(x) promoters(x, up = up, down = down))

# remove out of bound ranges
  domTSS_win.grl <- lapply(domTSS_win.grl, function(x) x[width(trim(x)) == win])

# check the length of the each range in the list
  lapply(domTSS_win.grl, length)

# randomize prior to sorting
  set.seed(10)
  random_row.l <- lapply(sapply(domTSS_win.grl, length), sample)

  for (i in 1:length(domTSS_win.grl)) {
    domTSS_win.grl[[i]] <- domTSS_win.grl[[i]][random_row.l[[i]]]  
  }

# sort index - by the distance and orientation of the shift
  sort_idx.l <- sapply(domTSS_win.grl, function(x) order(x$domTSS_dist,  decreasing = FALSE))
  for (i in 1:length(domTSS_win.grl)) {
    domTSS_win.grl[[i]] <- domTSS_win.grl[[i]][sort_idx.l[[i]]]  
  }

# extract sequence 
  domTSS_seq.l <- sapply(domTSS_win.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))


# plot dinucleotide plots arround nucleosome centers (ranges are centered on nucleosome dyads)
  pattern <- c("WW", "SS", "TATAWAWR")
  hm.l <- list()
  hm_smoothed.l <- list()
  scale_range <- list(c(0, 0.6), c(0, 0.8), c(0, 0.068))
  names(scale_range) <- pattern

for (i in 1:length(pattern)) {
  for(j in 1:length(domTSS_seq.l)) {
    hm.l = PatternHeatmap(domTSS_seq.l[[j]], pattern = pattern[i], coords = range, label = paste(pattern[i], ", ", samples[j], sep = ""))
    hm_smoothed.l = smoothHeatmap(hm.l, sigma = c(3, 3), output.size=c(500, 500))
    scale(hm_smoothed.l) <- scale_range[[i]]
    
    pdf(paste("Figures/heatmaps/shifting_promoters/shift_vs_mESC/domTSS_sample_mESC_cent", pattern[[i]], samples[[j]], ".pdf", sep = "_"), height = 4, width = 4.5)
    plotHeatmapList(hm_smoothed.l, legend = TRUE, cex.label = 0.5)
    dev.off()
  }}

## Plotting H3K4me3 patterns - MUU H3K4me3 ChIP example
  # calculate subtracted coverage 
    # load H3K4me3 object
      MII_H3K4me3_bam_1.gr <- readRDS("merged_replicates/intermediate_data/Liu_MII_H3K4me3_bam1_gr.RDS")
      MII_H3K4me3_bam_2.gr <- readRDS("merged_replicates/intermediate_data/Liu_MII_H3K4me3_bam2_gr.RDS")
      MII_H3K4me3_bam_3.gr <- readRDS("merged_replicates/intermediate_data/Liu_MII_H3K4me3_bam3_gr.RDS")
      MII_H3K4me3_bam_4.gr <- readRDS("merged_replicates/intermediate_data/Liu_MII_H3K4me3_bam4_gr.RDS")
  
    # calculate coverage separately for + and - strand
      MII_H3K4me3_plus.cov1 <- coverage(MII_H3K4me3_bam_1.gr[as.character(strand(MII_H3K4me3_bam_1.gr)) == "+", ])
      MII_H3K4me3_minus.cov1 <- coverage(MII_H3K4me3_bam_1.gr[as.character(strand(MII_H3K4me3_bam_1.gr)) == "-", ])
      
      # calculate coverage separately for + and - strand
      MII_H3K4me3_plus.cov2 <- coverage(MII_H3K4me3_bam_2.gr[as.character(strand(MII_H3K4me3_bam_2.gr)) == "+", ])
      MII_H3K4me3_minus.cov2 <- coverage(MII_H3K4me3_bam_2.gr[as.character(strand(MII_H3K4me3_bam_2.gr)) == "-", ])  
      
      # calculate coverage separately for + and - strand
      MII_H3K4me3_plus.cov3 <- coverage(MII_H3K4me3_bam_3.gr[as.character(strand(MII_H3K4me3_bam_3.gr)) == "+", ])
      MII_H3K4me3_minus.cov3 <- coverage(MII_H3K4me3_bam_3.gr[as.character(strand(MII_H3K4me3_bam_3.gr)) == "-", ])  
      
      # calculate coverage separately for + and - strand
      MII_H3K4me3_plus.cov4 <- coverage(MII_H3K4me3_bam_4.gr[as.character(strand(MII_H3K4me3_bam_4.gr)) == "+", ])
      MII_H3K4me3_minus.cov4 <- coverage(MII_H3K4me3_bam_4.gr[as.character(strand(MII_H3K4me3_bam_4.gr)) == "-", ])  
      
  
  
    # calculated subtracted coverage (plus-minus strand)
      MII_H3K4me3_subtracted.cov1 <- MII_H3K4me3_plus.cov1 - MII_H3K4me3_minus.cov1
      MII_H3K4me3_subtracted.cov2 <- MII_H3K4me3_plus.cov2 - MII_H3K4me3_minus.cov2
      MII_H3K4me3_subtracted.cov3 <- MII_H3K4me3_plus.cov3 - MII_H3K4me3_minus.cov3
      MII_H3K4me3_subtracted.cov4 <- MII_H3K4me3_plus.cov4 - MII_H3K4me3_minus.cov4
    
      MII_H3K4me3_subtracted.cov <- MII_H3K4me3_subtracted.cov1 + MII_H3K4me3_subtracted.cov2 + MII_H3K4me3_subtracted.cov3 + MII_H3K4me3_subtracted.cov4
    
  
      samples <- names(domTSS_oocyte_som6.grl)  
    
    # heatmaps plotting
    # create windows for plotting
      up <- 500
      down <- 500
      range <- c(-up, down)
      win <- up + down
      
      # centre on domTSS 
        domTSS_win.grl <- lapply(domTSS_oocyte_som6.grl, function(x) {
          gr <- promoters(x, up = up, down = down)
        })
    
      # remove out of bound ranges
        domTSS_win.grl <- lapply(domTSS_win.grl, function(x) x[width(trim(x)) == win])
      
      # check the length of the each range in the list
        sapply(domTSS_win.grl, length)
      
      # randomize prior to sorting
        set.seed(10)
        random_row.l <- lapply(lapply(domTSS_win.grl, length), sample)
        
        for (i in 1:length(domTSS_win.grl)) {
          domTSS_win.grl[[i]] <- domTSS_win.grl[[i]][random_row.l[[i]]]  
        }
      
      # sort index - no sequences are trimmed out so I can use the .gr object for sorting
        sort_idx.l <- sapply(domTSS_win.grl, function(x) order(x$interquantile_width,  decreasing = FALSE))
        for (i in 1:length(domTSS_win.grl)) {
          domTSS_win.grl[[i]] <- domTSS_win.grl[[i]][sort_idx.l[[i]]]  
        }
      
        samples <- names(domTSS_win.grl)
      
      # coverage with H3K4me3 MII oocyte
        hm.l <- list()
        for (i in 1:length(domTSS_win.grl)) {
          hm.l[[i]] <- CoverageHeatmap(domTSS_win.grl[[i]],
                                       track =  MII_H3K4me3_subtracted.cov,
                                       coords = range,
                                       label = paste("H3K4me3 MII, ", samples[[i]], sep = ""))
          image(hm.l[[i]])[as.character(strand(domTSS_win.grl[[i]])) == "-", ] <- -image(hm.l[[i]])[as.character(strand(domTSS_win.grl[[i]])) == "-", ]
          
          scale(hm.l[[i]]) <- quantile(hm.l[[i]]@image, c(0.01, 0.99))
        }
        
        for(i in 1:length(hm.l)) {
          pdf(paste("Figures/heatmaps/H3K4me3/oocyte_spec_domTSS_H3K4me3_MII_Liu_subcov_all_", samples[[i]], ".pdf", sep = "_"), height = 5.5, width = 4.5)
          plotHeatmapList(hm.l[[i]],
                          cex.label=0.5,
                          color = c("blue", "white", "red"),
                          legend = TRUE)
          
          dev.off()
        }
