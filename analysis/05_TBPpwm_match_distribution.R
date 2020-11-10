#/*==========================================================================#*/
#' ## Figure 2G, 3E
#/*==========================================================================#*/

# code is the same for all graphs - different GRanges object can be used - different set of sequences and centring point
# load libraries
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(seqPattern)
  library(matrixStats)
  library(tidyr)
  library(ggplot2)

# load domTSS annotated GRanges - centred on sample domTSS
  domTSS_PGCs_oocyte_embryo_anno.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/domTSS_PGCs_oocyte_embryo_anno_grl.RDS")

#--- selecting W-box or TATA-box regions: -35 to -20 bp upstream of dominant TSS ---#
# create windows for plotting
  up <- 15
  down <- 0
  range <- c(-up, down)
  win <- up + down

  domTSS_win.grl <- lapply(domTSS_PGCs_oocyte_embryo_anno.grl, function(x) promoters(x, up = up, down = down))
  domTSS_win.grl <- lapply(domTSS_win.grl, function(x) shift(x, shift = -20*as.integer(ifelse(strand(x) =="-", -1, 1))))

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
  rowmax.l <- lapply(tbp_scores.l, rowMaxs)

# convert each list to a dataframe
  rowmax.dfl <- lapply(1:length(rowmax.l), function(x) {
    df <- data.frame("value" = rowmax.l[[x]], "sample" = names(rowmax.l)[x])
  })

# change to ggplot compatible
  rowmax.gg <- as.data.frame(do.call(rbind, rowmax.dfl))

# change levels
  rowmax.gg$sample <- factor(rowmax.gg$sample, levels = names(domTSS_PGCs_oocyte_embryo_anno.grl))

# overlay distributions for plotting
## select E16_5F and E14 mESC
  library(dplyr)
  library(MASS)

  rowmax.selectA.gg <- as_tibble(rowmax.gg) %>% filter(sample  == "E14_mESC" | sample  == "E16_5F") 

# drop unused levels
  rowmax.selectA.gg <- droplevels(rowmax.selectA.gg)

# plot distribution of domTSS distances 
  col <- col_scheme[c("mESC-E14", "E16.5F")]

# calculate medians for each distribution
  mESC_med <- median(rowmax.selectA.gg %>% filter(sample =="E14_mESC") %>% dplyr::select(value) %>% unlist)
  E16_5F_med <- median(rowmax.selectA.gg %>% filter(sample =="E16_5F") %>% dplyr::select(value) %>% unlist)

# plot
  plotA <- ggplot(rowmax.selectA.gg, aes(x = value,  fill = sample, alpha = 0.85)) +
    geom_histogram(aes(x =rowmax.selectA.gg$value , y = ..count../sum(..count..)*100), col = "black", lwd = 0.5, binwidth = 2, position = "identity", alpha = 0.4) +
    scale_fill_manual("sample", values = c("firebrick1", "#3FBC73FF")) +
    theme_light()  +
    theme(text = element_text(size = 20, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          legend.position = "top") +
    labs(y = "%", x = "TBP pwm percentile match") + 
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
    scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 1)) +
    geom_vline(xintercept = mESC_med, lty = "dashed", lwd = 1, col = col[1]) +
    geom_vline(xintercept = E16_5F_med, lty = "dashed", lwd = 1, col = col[2])
