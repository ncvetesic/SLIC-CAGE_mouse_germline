#/*==========================================================================#*/
#' ## Figure 2G, 3E
#/*==========================================================================#*/
# example code for seqologo generation - any domTSS Granges can be used, here zebrafish domTSS Granges object is used (as in Figure 2h)

#--- selecting W-box regions - centered on maternal TSS ---#
# create windows for plotting
  up <- 35
  down <- 5
  range <- c(-up, down)
  win <- up + down

  domTSS_win.gr <- promoters(domTSS_egg_p6_shift_fish_egg_centred.grl[[1]], up = up, down = down)

# extract sequence 
  domTSS_seq <- getSeq(BSgenome.Drerio.UCSC.danRer7, domTSS_win.gr)

# create consensus out of all sequences - seqlogos
library(ggseqlogo)
  pdf("zebrafish/seqlogo_shifting_egg_vs_prim6_egg_centred.pdf", height = 6, width = 6)  
    print(ggseqlogo(as.character(domTSS_seq)) +
            scale_y_continuous(limits = c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)))
  dev.off()

#--- Zoomed W-box - selecting W-box regions - centered on maternal TSS ---#
# create windows for plotting
  up <- 0
  down <- 7
  range <- c(-up, down)
  win <- up + down

  domTSS_win.gr <- promoters(domTSS_egg_p6_shift_fish_egg_centred.grl[[1]], up = up, down = down)
  domTSS_win.gr <- shift(domTSS_win.gr, shift = -30*as.integer(ifelse(strand(domTSS_win.gr) =="-", -1, 1)))

# extract sequence 
  domTSS_seq <- getSeq(BSgenome.Drerio.UCSC.danRer7, domTSS_win.gr)

# create consensus out of all sequences - seqlogos
  library(ggseqlogo)
  pdf("zebrafish/seqlogo_shifting_egg_vs_prim6_egg_centred_wbox_zoom.pdf", height = 6, width = 6)  
    print(ggseqlogo(as.character(domTSS_seq)))
  dev.off()