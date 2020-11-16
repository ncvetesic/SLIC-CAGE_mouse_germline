#/*==========================================================================#*/
#' ## Extended Data Figure 7A-D
#/*==========================================================================#*/
# Tetranucleotide analysis in -30 bp regions - shifting promoters

# load libraries
library(BSgenome.Mmusculus.UCSC.mm10)
library(seqPattern)
library(matrixStats)
library(tidyr)
library(ggplot2)
library(CAGEr)
library(ggseqlogo)
data(TBPpwm) #TBP pwm is from CAGEr

# load mouse shifting promoters centred on egg domTSS (list)
  domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/domTSS_PGC_oocyte_embryo_vs_mESC_shift_shift_X_grl.RDS") # y is MESC E14

#--- selecting W-box or TATA-box regions: -34 to -23 bp upstream of dominant TSS ---#
# create windows for plotting
  up <- 0
  down <- 12
  range <- c(-up, down)
  win <- up + down

  domTSS_win.grl <- lapply(domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl, function(x) promoters(x, up = up, down = down))
  domTSS_win.grl <- lapply(domTSS_win.grl, function(x) shift(x, shift = -34*as.integer(ifelse(strand(x) =="-", -1, 1))))

# remove out of bound ranges
  idx.l <- lapply(domTSS_win.grl, function(x) width(trim(x)) == win)
  for(i in 1:length(domTSS_win.grl)) {
    domTSS_win.grl[[i]] <- domTSS_win.grl[[i]][idx.l[[i]]]
  }

# extract sequence 
  domTSS_seq.l <- lapply(domTSS_win.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))

# scan with TBP motif
  tbp_scores.l <- lapply(domTSS_seq.l, function(x) motifScanScores(x, motifPWM = TBPpwm))

# select highest score per sequence
#rowmax.l <- lapply(tbp_scores.l, rowMaxs)

# select position with highest TBPpwm match per sequence
  tbp_positions.l <- lapply(1:length(tbp_scores.l), function(x) {
    df <- data.frame(sequence = 1:nrow(tbp_scores.l[[x]]), position = max.col(tbp_scores.l[[x]]))
  })

  tbp.grl <- list()
  for(i in 1:length(tbp_positions.l)) {
    tbp.grl[[i]] <- GRanges()
    for(j in 1:nrow(tbp_positions.l[[i]])) {
      tbp.grl[[i]] <- c(tbp.grl[[i]], GRanges(seqnames = seqnames(domTSS_win.grl[[i]][tbp_positions.l[[i]][j, "sequence"]]),
                                              ranges = IRanges(start = start(domTSS_win.grl[[i]][tbp_positions.l[[i]][j, "sequence"]]) + tbp_positions.l[[i]][j, "position"],
                                                               end = start(domTSS_win.grl[[i]][tbp_positions.l[[i]][j, "sequence"]]) + tbp_positions.l[[i]][j, "position"] + 7),
                                              strand = strand(domTSS_win.grl[[i]][tbp_positions.l[[i]][j, "sequence"]])))
    }}


  names(tbp.grl) <- names(tbp_positions.l)

# extract sequence 
  tbp_seq.l <- lapply(tbp.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))
  
  names(tbp_seq.l) <- names(domTSS_win.grl)

# plot sequence logos of highest scoring sequences
  library(ggseqlogo)
  for(i in 1:length(tbp_seq.l)) {
    pdf(paste0("../all_replicates/merged_replicates/results/wbox_seqlogos/TATA_like_", names(tbp_seq.l)[i], ".pdf"), height = 6, width = 8)  
    print(ggseqlogo(as.character(tbp_seq.l[[i]])) + 
            scale_y_continuous(limits = c(0, 0.2), breaks=c(0, 0.05, 0.1, 0.15, 0.2)))
    dev.off()  
  }


# tetramer count in TATA regions
# split into characters
  tbp_seq_char.l <- lapply(tbp_seq.l, function(x) {
    x <- as.character(x)
    x <- strsplit(x, "")
    return(x) 
  })
  names(tbp_seq_char.l) <- names(tbp_seq.l)

# select only oocyte samples
  tbp_seq_char_sel.l <- tbp_seq_char.l[c("PN6", "oocyte", "S4_cell")]

  tetramer_res.l <- lapply(tbp_seq_char_sel.l, tetramer_count)

# convert to dataframe
  samples <- names(tetramer_res.l)
  tetramers.df <- list()

  for(i in 1:length(tetramer_res.l)) {
    tetramers.df[[i]] <- as.data.frame(tetramer_res.l[[i]])
    tetramers.df[[i]]$sample <- samples[i]
  }
  tetramers.gg <- do.call(rbind, tetramers.df)
  
  tetramers_top5.df <- lapply(tetramers.df, function(x){
    return(x[1:10, ])
  })

  tetramers_top5.gg <- do.call(rbind, tetramers_top5.df)

# plot tetramer frequency
  tetramers_top5.gg$sample <- factor(tetramers_top5.gg$sample, levels = samples)
  tetramers_top5.gg$tetramers <- droplevels(tetramers_top5.gg$tetramers)
  tetramers_top5.gg$tetramers <- factor(tetramers_top5.gg$tetramers, levels = rev(sort(levels(tetramers_top5.gg$tetramers))))  

  rownames(tetramers_top5.gg) <- NULL

  library(ggplot2)
  library(viridis)
  col = rep("gray66", times = 20)

  p <- ggplot(data = tetramers_top5.gg, aes(x = tetramers,
                                            y = Freq, fill = tetramers)) + 
    scale_fill_manual(values = col) +
    geom_bar(stat = "identity", position = position_dodge(), colour = "black", size = 0.25) + 
    coord_flip() + 
    scale_y_continuous(name = "%", limits = c(0, 2), seq(0, 2, by = 0.5)) +
    xlab("TATA-like regions") + ggtitle("Tetramers in TATA-like regions") + 
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12))
  
  pdf("Figures/TATA_tetramers_perc_domTSS_egg.pdf", height = 4, width = 8)
  p + facet_wrap(~sample, ncol = 3, scales = "free")
  dev.off()

### --- repeat with centring on somatic dominant TSS --- ###
# load mouse shifting promoters centred on mESC domTSS (list)
  domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/domTSS_PGC_oocyte_embryo_vs_mESC_shift_grl.RDS") # y is MESC E14

#--- selecting W-box or TATA-box regions: -34 to -23 bp upstream of dominant TSS ---#
# create windows for plotting
  up <- 0
  down <- 12
  range <- c(-up, down)
  win <- up + down

  domTSS_win.grl <- lapply(domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl, function(x) promoters(x, up = up, down = down))
  domTSS_win.grl <- lapply(domTSS_win.grl, function(x) shift(x, shift = -34*as.integer(ifelse(strand(x) =="-", -1, 1))))

# remove out of bound ranges
  idx.l <- lapply(domTSS_win.grl, function(x) width(trim(x)) == win)
  for(i in 1:length(domTSS_win.grl)) {
    domTSS_win.grl[[i]] <- domTSS_win.grl[[i]][idx.l[[i]]]
  }

# extract sequence 
  domTSS_seq.l <- lapply(domTSS_win.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))

# scan with TBP motif
  tbp_scores.l <- lapply(domTSS_seq.l, function(x) motifScanScores(x, motifPWM = TBPpwm))

# select position with highest TBPpwm match per sequence
  tbp_positions.l <- lapply(1:length(tbp_scores.l), function(x) {
    df <- data.frame(sequence = 1:nrow(tbp_scores.l[[x]]), position = max.col(tbp_scores.l[[x]]))
  })

  tbp.grl <- list()
  for(i in 1:length(tbp_positions.l)) {
    tbp.grl[[i]] <- GRanges()
    for(j in 1:nrow(tbp_positions.l[[i]])) {
      tbp.grl[[i]] <- c(tbp.grl[[i]], GRanges(seqnames = seqnames(domTSS_win.grl[[i]][tbp_positions.l[[i]][j, "sequence"]]),
                                              ranges = IRanges(start = start(domTSS_win.grl[[i]][tbp_positions.l[[i]][j, "sequence"]]) + tbp_positions.l[[i]][j, "position"],
                                                               end = start(domTSS_win.grl[[i]][tbp_positions.l[[i]][j, "sequence"]]) + tbp_positions.l[[i]][j, "position"] + 7),
                                              strand = strand(domTSS_win.grl[[i]][tbp_positions.l[[i]][j, "sequence"]])))
    }}


  names(tbp.grl) <- names(tbp_positions.l)

# extract sequence 
  tbp_seq.l <- lapply(tbp.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))

  names(tbp_seq.l) <- names(domTSS_win.grl)

  library(ggseqlogo)
  for(i in 1:length(tbp_seq.l)) {
    pdf(paste0("../all_replicates/merged_replicates/results/wbox_seqlogos/TATA_like_mESC_centred_", names(tbp_seq.l)[i], ".pdf"), height = 6, width = 8)  
    print(ggseqlogo(as.character(tbp_seq.l[[i]])) + 
            scale_y_continuous(limits = c(0, 0.2), breaks=c(0, 0.05, 0.1, 0.15, 0.2)))
    dev.off()  
  }


# tetramer count in TATA regions
  # split into characters
  tbp_seq_char.l <- lapply(tbp_seq.l, function(x) {
    x <- as.character(x)
    x <- strsplit(x, "")
    return(x) 
  })
  names(tbp_seq_char.l) <- names(tbp_seq.l)

# select only oocyte samples
  tbp_seq_char_sel.l <- tbp_seq_char.l[c("PN6", "oocyte", "S4_cell")]

  tetramer_res.l <- lapply(tbp_seq_char_sel.l, tetramer_count)

# convert to dataframe
  samples <- names(tetramer_res.l)

  tetramers.df <- list()

  for(i in 1:length(tetramer_res.l)) {
    tetramers.df[[i]] <- as.data.frame(tetramer_res.l[[i]])
    tetramers.df[[i]]$sample <- samples[i]
  }
  tetramers.gg <- do.call(rbind, tetramers.df)

  tetramers_top5.df <- lapply(tetramers.df, function(x){
    return(x[1:10, ])
  })

  tetramers_top5.gg <- do.call(rbind, tetramers_top5.df)

# plot tetramer frequency
  tetramers_top5.gg$sample <- factor(tetramers_top5.gg$sample, levels = samples)
  tetramers_top5.gg$tetramers <- droplevels(tetramers_top5.gg$tetramers)
  tetramers_top5.gg$tetramers <- factor(tetramers_top5.gg$tetramers, levels = rev(sort(levels(tetramers_top5.gg$tetramers))))  

  rownames(tetramers_top5.gg) <- NULL

  library(ggplot2)
  col = rep("gray66", times = 20)

  p <- ggplot(data = tetramers_top5.gg, aes(x = tetramers,
                                            y = Freq, fill = tetramers)) + 
    scale_fill_manual(values = col) +
    geom_bar(stat = "identity", position = position_dodge(), colour = "black", size = 0.25) + 
    coord_flip() + 
    scale_y_continuous(name = "%", limits = c(0, 2), seq(0, 2, by = 0.5)) +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12))
  
    pdf("Figures/TATA_tetramers_perc_domTSS_som.pdf", height = 4, width = 8)
      p + facet_wrap(~sample, ncol = 3, scales = "free")
    dev.off()

