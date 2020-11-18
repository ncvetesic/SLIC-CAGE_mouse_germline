#/*==========================================================================#*/
#' ## Extended Data Figure 13A-D
#/*==========================================================================#*/
# Generation of nucleosome pre-positioning pwm and analyses
## PWM for nucleosome positioning and all metaplots and heatmaps

library(BSgenome.Mmusculus.UCSC.mm10)
library(magrittr)
library(ggseqlogo)

# load domTSS object
  domTSS_PGCs_oocyte_embryo_anno.grl <- readRDS(file ="../all_reps_PN6_2cell/intermediate_data/domTSS_PGCs_oocyte_embryo_anno_grl.RDS")
  samples <- names(domTSS_PGCs_oocyte_embryo_anno.grl)

# WW extraction on [6 to 50 bp] and no smoothing, binsize = 1
# create windows centered on domTSS 
  up <- 0
  down <- 45
  range <- c(-up, down)
  win <- up + down

# select only mESC
  domTSS_E14_win.gr <- domTSS_PGCs_oocyte_embryo_anno.grl[["E14_mESC"]]
  domTSS_E14_win.gr <- promoters(domTSS_E14_win.gr, up = up, down = down)

# exclude E16_5F shifting promoters
  PGC_oocyte_embryo_vs_mESC_shift_anno.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/PGC_oocyte_embryo_vs_mESC_shift_anno_grl.RDS")

# find overlaps
  overlaps <- findOverlaps(PGC_oocyte_embryo_vs_mESC_shift_anno.grl[["E16_5F"]], domTSS_E14_win.gr)

# exclude overlaps
  domTSS_E14_win.gr <- domTSS_E14_win.gr[-subjectHits(overlaps)]
  domTSS_mESC.gr <- domTSS_PGCs_oocyte_embryo_anno.grl[[1]]
  domTSS_mESC.gr <- domTSS_mESC.gr[-subjectHits(overlaps)]

  domTSS_E14_win_shift.gr <- shift(domTSS_E14_win.gr, shift = 5*as.integer(ifelse(as.character(strand(domTSS_E14_win.gr)) == "-", -1, 1)))

# remove out of bound ranges
  domTSS_E14_win_shift.gr <- domTSS_E14_win_shift.gr[width(trim(domTSS_E14_win_shift.gr)) == win]

# separate promoters into broad and sharp
  broad_idx <- domTSS_E14_win_shift.gr$interquantile_width > 9
  sharp_idx <- domTSS_E14_win_shift.gr$interquantile_width <= 9

  domTSS_E14_broad_win.gr <- domTSS_E14_win_shift.gr[broad_idx]
  domTSS_E14_sharp_win.gr <- domTSS_E14_win_shift.gr[sharp_idx]

  domTSS_E14_broad.gr <- domTSS_mESC.gr[broad_idx]
  domTSS_E14_sharp.gr <- domTSS_mESC.gr[sharp_idx]

# select 10k broad promoters to create pwm
  set.seed(7)
  select <- sample(length(domTSS_E14_broad_win.gr), 8000)
  domTSS_E14_broad_sel_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, domTSS_E14_broad_win.gr[select])

# create nucleosome pwm
  library(Biostrings)
  nucleosome_pwm <- consensusMatrix(domTSS_E14_broad_sel_seq)
  nucleosome_pwm <- toPWM(nucleosome_pwm)

# save sequence composition of nucleosom pwm
  pdf("Figures/heatmaps/nucl_pwm_seqlogo.pdf", height = 3, width = 6)
    print(ggseqlogo(as.character(domTSS_E14_broad_sel_seq)))
  dev.off()

# save nucleosome_pwm as object
  saveRDS(nucleosome_pwm, "intermediate_data/nucleosome_pwm.RDS")

# use nucleosome pwm to scan leftover sequences for testing
# create windows for plotting
  up <- 1000
  down <- 1000
  range <- c(-up, down)
  win <- up + down

# centre on domTSS 
  domTSS_win.grl <- lapply(domTSS_PGCs_oocyte_embryo_anno.grl, function(x) {
    gr <- promoters(x, up = up, down = down)
  })

# remove out of bound ranges
  domTSS_win.grl <- lapply(domTSS_win.grl, function(x) x[width(trim(x)) == win])
  domTSS_win.grl <- domTSS_win.grl[[1]]

# sort index - no sequences are trimmed out so I can use the .gr object for sorting
  sort_idx.l <- order(domTSS_win.grl$interquantile_width,  decreasing = FALSE)
  domTSS_win.grl <- domTSS_win.grl[sort_idx.l]  

  domTSS_win_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, domTSS_win.grl)

# scan the sequence with nucleosome pwm
  library(seqPattern)
  nucl_scores <- motifScanScores(domTSS_win_seq, motifPWM = nucleosome_pwm)  

# use heatmaps for plotting
  hm.l <- new("Heatmap",
              image = nucl_scores,
              coords = as.integer(range),
              nseq = nrow(nucl_scores),
              metadata = list())
  
  hm_smoothed.l<- smoothHeatmap(hm.l, sigma = c(3, 3), output.size=c(500, 500))
  scale(hm_smoothed.l) <- c(40, 60)

# plot heatmaps
  pdf("Figures/heatmaps/nucl_pwm_mESC_sharp_broad.pdf", height = 5.5, width = 4.5)
  plotHeatmapList(hm_smoothed.l,
                  legend = T,
                  legend.pos = "l",
                  color = "Spectral",
                  cex.label = 0.5)
  dev.off()


# plot metaplots - for test promoters, not used in pwm creation
  domTSS_broad_test.gr <- domTSS_E14_broad.gr[-select]

# select 825 random promoters
  select2 <- sample(length(domTSS_broad_test.gr), size = length(domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl[["E16_5F"]]))
  domTSS_broad_test.gr <- domTSS_broad_test.gr[select2]

# use nucleosome pwm to scan leftover sequences for testing
# create windows for plotting
  up <- 150
  down <- 500
  range <- c(-up, down)
  win <- up + down

# centre on domTSS 
  domTSS_win.grl <- promoters(domTSS_broad_test.gr, up = up, down = down)

# remove out of bound ranges
  domTSS_win.grl <- domTSS_win.grl[width(trim(domTSS_win.grl)) == win]

# sort index - no sequences are trimmed out so I can use the .gr object for sorting
sort_idx.l <- order(domTSS_win.grl$interquantile_width,  decreasing = FALSE)
domTSS_win.grl <- domTSS_win.grl[sort_idx.l]  

# extract sequence
  domTSS_broad_test_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, domTSS_win.grl)

# scan the sequence with nucleosome pwm
  library(seqPattern)
  nucl_scores <- motifScanScores(domTSS_broad_test_seq, motifPWM = nucleosome_pwm)  

# produce heatmaps prior to plotting (no smoothing)
  hm_broad.l <- list()

# use heatmaps for plotting
  hm_broad.l <- new("Heatmap",
                    image = nucl_scores,
                    coords = as.integer(range),
                    nseq = nrow(nucl_scores),
                    metadata = list())
  library(DescTools) 
  matrix_sum <- colSums(Winsorize(hm_broad.l@image))/nrow(hm_broad.l@image)

# create a list of dataframes (with sharp and broad column)
  all_df_test <- data.frame("x_coord" = -up:(down-ncol(nucleosome_pwm)), 
                            "broad_occurrence" = matrix_sum,
                            "sample" = "broad_test")


# plot metaplots using ggplot2
  library(ggplot2)
  p <- ggplot(all_df_test) +
    geom_line(lty = 1, aes(x = x_coord, y = broad_occurrence), size = 0.6) +
    scale_color_manual(values = "black", 
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
    scale_y_continuous(name = "Relative score", limits = c(50, 60), breaks = seq(50, 60, by = 2)) +
    scale_x_continuous(name = "Distance to dominant TSS (bp)", limits = c(-150, 350), seq(-150, 350, by = 50))

    pdf("Figures/heatmaps/nucl_pwm_test_broad_meta.pdf", height = 4, width = 6)
      print(p)
    dev.off()


# test E16.5F shifting promoters (X is sample centered, no X is mESC centred)
  domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/domTSS_PGC_oocyte_embryo_vs_mESC_shift_grl.RDS")
  domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/domTSS_PGC_oocyte_embryo_vs_mESC_shift_shift_X_grl.RDS")

# select E16.5 stage
  domTSS_E16_5F_shift_mESC.gr <- domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl[["E16_5F"]]
  domTSS_E16_5F_shift_E16_5F.gr <- domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl[["E16_5F"]]

# create windows for plotting
  up <- 150
  down <- 500
  range <- c(-up, down)
  win <- up + down

# centre on domTSS 
  domTSS_win.grl <- promoters(domTSS_E16_5F_shift_E16_5F.gr, up = up, down = down)

# remove out of bound ranges
  domTSS_win.grl <- domTSS_win.grl[width(trim(domTSS_win.grl)) == win]

# extract sequence
  domTSS_broad_test_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, domTSS_win.grl)

# scan the sequence with nucleosome pwm
  library(seqPattern)
  nucl_scores <- motifScanScores(domTSS_broad_test_seq, motifPWM = nucleosome_pwm)  

# produce heatmaps prior to plotting (no smoothing)
  hm_broad.l <- list()

# use heatmaps for plotting
  hm_broad.l <- new("Heatmap",
                    image = nucl_scores,
                    coords = as.integer(range),
                    nseq = nrow(nucl_scores),
                    metadata = list())


  matrix_sum <- colSums(Winsorize(hm_broad.l@image))/nrow(hm_broad.l@image)


# create a list of dataframes (with sharp and broad column)
  all_df_E16_5F <- data.frame("x_coord" = -up:(down-ncol(nucleosome_pwm)), 
                              "broad_occurrence" = matrix_sum,
                              "sample" = "E16.5_domTSS")

# create windows for plotting
  up <- 150
  down <- 500
  range <- c(-up, down)
  win <- up + down

# centre on domTSS 
  domTSS_win.grl <- promoters(domTSS_E16_5F_shift_mESC.gr, up = up, down = down)

# remove out of bound ranges
  domTSS_win.grl <- domTSS_win.grl[width(trim(domTSS_win.grl)) == win]

# extract sequence
  domTSS_broad_test_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, domTSS_win.grl)

# scan the sequence with nucleosome pwm
  library(seqPattern)
  nucl_scores <- motifScanScores(domTSS_broad_test_seq, motifPWM = nucleosome_pwm)  

# produce heatmaps prior to plotting (no smoothing)
  hm_broad.l <- list()
  
  # use heatmaps for plotting
  hm_broad.l <- new("Heatmap",
                    image = nucl_scores,
                    coords = as.integer(range),
                    nseq = nrow(nucl_scores),
                    metadata = list())


  matrix_sum <- colSums(Winsorize(hm_broad.l@image))/nrow(hm_broad.l@image)

# create a list of dataframes (with sharp and broad column)
  all_df_mESC <- data.frame("x_coord" = -up:(down-ncol(nucleosome_pwm)), 
                            "broad_occurrence" = matrix_sum,
                            "sample" = "mESC_domTSS")
  
  # combine 3 dataframes (centring on domTSS from mESC or E16_5F and broad test)
    all_df <- rbind(all_df_test, all_df_mESC, all_df_E16_5F)

# change levels
  all_df$sample <- factor(all_df$sample, levels = c("broad_test", "mESC_domTSS", "E16.5_domTSS"))

# plot metaplots using ggplot2 - all distances included
  library(ggplot2)
  p <- ggplot(all_df) +
    geom_line(lty = 1, aes(x = x_coord, y = broad_occurrence, colour = sample), size = 0.6) +
    scale_color_manual(values = c("black", "firebrick1", "#3FBC73FF"), 
                       name = NULL) +
    theme_light()   +
    theme(text = element_text(size = 14, colour = "black"),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          legend.position = "top",
          legend.title = element_blank()) +
    guides(colour = guide_legend(reverse=F)) +
    geom_vline(xintercept = 0, lty = "dashed", colour = "gray48", size = 0.6) +
    scale_y_continuous(name = "Relative score", limits = c(56, 64), breaks = seq(56, 64, by =  2)) +
    scale_x_continuous(name = "Distance to dominant TSS (bp)", limits = c(-150, 200), seq(-150, 200, by = 25))
  
  
  
  pdf("Figures/heatmaps/nucl_pwm_all.pdf", height = 4, width = 5)
    print(p)
  dev.off()


# use E16_5F shifting promoters where there is a physical shift - select certain distance between domTSS and separate metaplots
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(seqPattern)

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

# filter to include shifting promoters that moved by at least 40 bp - other plots are recapitulated same way just changed distance
  domTSS_overlap.l <- lapply(domTSS_overlap.l, function(x) {
    sample <- x[x$domTSS_dist <= -40]
    return(sample)
  })

# select only E16_5F
  domTSS_overlap_E16_5F_X.gr <- domTSS_overlap.l[["E16_5F"]]

# overlap shifting promoters and tc to get iq-width/annotate with mcols
  domTSS_overlap.l <- list()
  for(i in 1:length(domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl)) {
    overlap <- findOverlapPairs(domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl[[i]], tc_PGCs_oocyte_embryo_anno_nomESC.grl[[i]])
    domTSS_overlap.l[i] <- overlap@first
    mcols(domTSS_overlap.l[[i]]) <- cbind(mcols(overlap@first), mcols(overlap@second))
  }
  names(domTSS_overlap.l) <- names(domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl)

# filter to include shifting promoters that moved by at least one bp
  domTSS_overlap.l <- lapply(domTSS_overlap.l, function(x) {
    sample <- x[x$domTSS_dist <= -40]
    return(sample)
  })

  domTSS_overlap_E16_5F.gr <- domTSS_overlap.l[["E16_5F"]]

# -----plotting

# create windows for plotting
  up <- 150
  down <- 500
  range <- c(-up, down)
  win <- up + down

# centre on domTSS 
  domTSS_win.grl <- promoters(domTSS_overlap_E16_5F_X.gr, up = up, down = down)

# remove out of bound ranges
  domTSS_win.grl <- domTSS_win.grl[width(trim(domTSS_win.grl)) == win]

# extract sequence
  domTSS_broad_test_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, domTSS_win.grl)
  
# scan the sequence with nucleosome pwm
  library(seqPattern)
  nucl_scores <- motifScanScores(domTSS_broad_test_seq, motifPWM = nucleosome_pwm)  

# produce heatmaps prior to plotting (no smoothing)
  hm_broad.l <- list()

# use heatmaps for plotting
  hm_broad.l <- new("Heatmap",
                    image = nucl_scores,
                    coords = as.integer(range),
                    nseq = nrow(nucl_scores),
                    metadata = list())


  matrix_sum <- colSums(Winsorize(hm_broad.l@image, probs = c(0.05, 0.95)))/nrow(hm_broad.l@image)

# create a list of dataframes (with sharp and broad column)
  all_df_E16_5F <- data.frame("x_coord" = -up:(down-ncol(nucleosome_pwm)), 
                              "broad_occurrence" = matrix_sum,
                              "sample" = "E16.5_domTSS")

# create windows for plotting
  up <- 150
  down <- 500
  range <- c(-up, down)
  win <- up + down

# centre on domTSS 
  domTSS_win.grl <- promoters(domTSS_overlap_E16_5F.gr, up = up, down = down)

# remove out of bound ranges
  domTSS_win.grl <- domTSS_win.grl[width(trim(domTSS_win.grl)) == win]

# extract sequence
  domTSS_broad_test_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, domTSS_win.grl)

# scan the sequence with nucleosome pwm
  library(seqPattern)
  nucl_scores <- motifScanScores(domTSS_broad_test_seq, motifPWM = nucleosome_pwm)  

# produce heatmaps prior to plotting (no smoothing)
  hm_broad.l <- list()

# use heatmaps for plotting
  hm_broad.l <- new("Heatmap",
                    image = nucl_scores,
                    coords = as.integer(range),
                    nseq = nrow(nucl_scores),
                    metadata = list())


  matrix_sum <- colSums(Winsorize(hm_broad.l@image, probs = c(0.05, 0.95)))/nrow(hm_broad.l@image)

# create a list of dataframes (with sharp and broad column)
  all_df_mESC <- data.frame("x_coord" = -up:(down-ncol(nucleosome_pwm)), 
                            "broad_occurrence" = matrix_sum,
                            "sample" = "mESC_domTSS")

# combine 3 dataframes (centring on domTSS from mESC or E16_5F and broad test)
  all_df <- rbind(all_df_test, all_df_mESC, all_df_E16_5F)

# change levels
  all_df$sample <- factor(all_df$sample, levels = c("broad_test", "mESC_domTSS", "E16.5_domTSS"))

# plot metaplots using ggplot2
  library(ggplot2)
  p <- ggplot(all_df) +
    geom_line(lty = 1, aes(x = x_coord, y = broad_occurrence, colour = sample), size = 0.6) +
    scale_color_manual(values = c("black", "firebrick1", "#3FBC73FF"), 
                       name = NULL) +
    theme_light()   +
    theme(text = element_text(size = 14, colour = "black"),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          legend.position = "top",
          legend.title = element_blank()) +
    guides(colour = guide_legend(reverse=F)) +
    geom_vline(xintercept = 0, lty = "dashed", colour = "gray48", size = 0.6) +
    scale_y_continuous(name = "Relative score", limits = c(56, 68), breaks = seq(56, 68, by =  2)) +
    scale_x_continuous(name = "Distance to dominant TSS (bp)", limits = c(-150, 200), seq(-150, 200, by = 25))



  pdf("Figures/heatmaps/nucl_pwm_all_40bp_upstream_shifted.pdf", height = 4, width = 5)
    print(p)
  dev.off()

# heatmap pwm nucleosome in shifting promoters
# create windows for plotting
  up <- 500
  down <- 500
  range <- c(-up, down)
  win <- up + down

# centre on mESC 
  domTSS_win.grl <- sapply(domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl, function(x) promoters(x, up = up, down = down))

# remove out of bound ranges
  domTSS_win.grl <- sapply(domTSS_win.grl, function(x) x[width(trim(x)) == win])

# check the length of the each range in the list
  sapply(domTSS_win.grl, length)

# randomize prior to sorting
  set.seed(10)
  random_row.l <- lapply(sapply(domTSS_win.grl, length), sample)

for (i in 1:length(domTSS_win.grl)) {
  domTSS_win.grl[[i]] <- domTSS_win.grl[[i]][random_row.l[[i]]]  
}

# sort index 
  sort_idx.l <- sapply(domTSS_win.grl, function(x) order(x$domTSS_dist,  decreasing = FALSE))
  for (i in 1:length(domTSS_win.grl)) {
    domTSS_win.grl[[i]] <- domTSS_win.grl[[i]][sort_idx.l[[i]]]  
  }

# extract sequence 
  domTSS_seq.l <- lapply(domTSS_win.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))

# select only 16_5F
  domTSS_seq.l <- domTSS_seq.l[["E16_5F"]]

# scan the sequence with nucleosome pwm
  library(seqPattern)
  nucl_scores <- motifScanScores(domTSS_seq.l, motifPWM = nucleosome_pwm)  

# produce heatmaps prior to plotting (no smoothing)
  hm_broad.l <- list()

# use heatmaps for plotting
  hm_broad.l <- new("Heatmap",
                    image = nucl_scores,
                    coords = as.integer(range),
                    nseq = nrow(nucl_scores),
                    metadata = list())

  hm_broad_mESC.l <- hm_broad.l

  hm_smoothed.l<- smoothHeatmap(hm_broad.l, sigma = c(3, 3), output.size=c(800, 1000))
  scale(hm_smoothed.l) <- c(40, 65)

# plot heatmaps
  pdf("Figures/heatmaps/shifting_promoters/shift_vs_mESC/E16_5_shifting_mESC_centred_nucl_pwm.pdf", height = 4, width = 4.5)
  plotHeatmapList(hm_smoothed.l,
                  legend = T,
                  legend.pos = "l",
                  color = "Spectral",
                  cex.label = 0.5)
  dev.off()


# centred on E16_5F
# create windows for plotting
  up <- 500
  down <- 500
  range <- c(-up, down)
  win <- up + down

# centre on E16_5F - as I dont see nucleosome shift signature (Mmusculus doesnt have it, I could see W-box shift)
  domTSS_win.grl <- sapply(domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl, function(x) promoters(x, up = up, down = down))

# remove out of bound ranges
  domTSS_win.grl <- sapply(domTSS_win.grl, function(x) x[width(trim(x)) == win])

# check the length of the each range in the list
  sapply(domTSS_win.grl, length)

# randomize prior to sorting
  set.seed(10)
  random_row.l <- lapply(sapply(domTSS_win.grl, length), sample)
  
  for (i in 1:length(domTSS_win.grl)) {
    domTSS_win.grl[[i]] <- domTSS_win.grl[[i]][random_row.l[[i]]]  
  }

# sort index
  sort_idx.l <- sapply(domTSS_win.grl, function(x) order(x$domTSS_dist,  decreasing = FALSE))
  for (i in 1:length(domTSS_win.grl)) {
    domTSS_win.grl[[i]] <- domTSS_win.grl[[i]][sort_idx.l[[i]]]  
  }

# extract sequence 
  domTSS_seq.l <- lapply(domTSS_win.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))

# select only 16_5F
  domTSS_seq.l <- domTSS_seq.l[["E16_5F"]]

# scan the sequence with nucleosome pwm
  library(seqPattern)
  nucl_scores <- motifScanScores(domTSS_seq.l, motifPWM = nucleosome_pwm)  

# produce heatmaps prior to plotting (no smoothing)
  hm_broad.l <- list()

# use heatmaps for plotting
  hm_broad.l <- new("Heatmap",
                    image = nucl_scores,
                    coords = as.integer(range),
                    nseq = nrow(nucl_scores),
                    metadata = list())

  hm_broad_E16_5F.l <- hm_broad.l

  hm_smoothed.l<- smoothHeatmap(hm_broad.l, sigma = c(3, 3), output.size=c(800, 1000))
  scale(hm_smoothed.l) <- c(40, 65)

# plot heatmaps
  pdf("Figures/heatmaps/shifting_promoters/shift_vs_mESC/E16_5_shifting_E16_5F_centred_nucl_pwm.pdf", height = 4, width = 4.5)
  plotHeatmapList(hm_smoothed.l,
                  legend = T,
                  legend.pos = "l",
                  color = "Spectral",
                  cex.label = 0.5)
  dev.off()