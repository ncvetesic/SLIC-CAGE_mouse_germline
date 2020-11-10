#/*==========================================================================#*/
#' ## Figure 2J
#/*==========================================================================#*/
# identifying W-box stretch length in mouse vs zebrafish
# prior to this analysis zebrafish CAGE data needs to be prepared/processed - see zebrafish_CAGE_shifts.R

library(seqinr)

# load mouse shifting promoters centred on egg domTSS (list)
  domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/domTSS_PGC_oocyte_embryo_vs_mESC_shift_shift_X_grl.RDS") # y is MESC E14

# load fish shifting promoters centred on egg domTSS
  domTSS_egg_p6_shift_fish_egg_centred.grl <- readRDS("zebrafish/domTSS_egg_p6_shift_fish_egg_centred_grl.RDS")


# --- zebrafish maternal centred counts --- #
# select 12 bp from -34 to -23
  up <- 0
  down <- 12
  range <- c(-up, down)
  win <- up + down

  domTSS_win.gr <- promoters(domTSS_egg_p6_shift_fish_egg_centred.grl[[1]], up = up, down = down)
  domTSS_win.gr <- shift(domTSS_win.gr, shift = -34*as.integer(ifelse(strand(domTSS_win.gr) =="-", -1, 1)))

# extract sequence 
  domTSS_seq <- getSeq(BSgenome.Drerio.UCSC.danRer7, domTSS_win.gr)
  domTSS.l <- as.list(as.character(domTSS_seq))
  domTSS_char.l <- lapply(domTSS.l, function(x) {
    s2c(x)
  })

  
#-----function to count WW stretches-----#
  ww_box_length <- function(seq, wordsize, start = 0, by = 1, freq = FALSE, alphabet = s2c("AT")){
    
    # set word sizes from 1 to set wordsize
    size_range <- 1:wordsize
    
    word_count.l <- lapply(size_range, function(x) {
      count_el <- count(seq = seq, wordsize = x, start = start, by = by,
                        freq = freq, alphabet = alphabet)
      
      # output is table - convert to named vector where elements are counts
      oligos <- names(count_el)
      count_el <- as.vector(count_el)
      names(count_el) <- oligos
      
      if(sum(count_el) >= 1) {
        max_word <- x
        return(x)
      } else {
        return(0)
      }
    })
    return(max(as.vector(unlist(word_count.l))))
  }
  
# count WW stretches
  word_size <- 12
  ww_box_stretch.l <- lapply(domTSS_char.l, function(x) {
    ww_box_length(x, word_size)
  })

  ww_box_stretch_fish <- as.vector(unlist(ww_box_stretch.l))

### --- mouse  maternal centred counts--- ###
# select 10 bp from -34 to -23
  up <- 0
  down <- 12
  range <- c(-up, down)
  win <- up + down

  domTSS_win.gr <- promoters(domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl[["PN6"]], up = up, down = down)
  domTSS_win.gr <- shift(domTSS_win.gr, shift = -34*as.integer(ifelse(strand(domTSS_win.gr) =="-", -1, 1)))
  
  # extract sequence 
  domTSS_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, domTSS_win.gr)
  domTSS.l <- as.list(as.character(domTSS_seq))
  domTSS_char.l <- lapply(domTSS.l, function(x) {
    s2c(x)
  })

# count WW stretches
  word_size <- 12
  ww_box_stretch.l <- lapply(domTSS_char.l, function(x) {
    ww_box_length(x, word_size)
  })

  ww_box_stretch_mouse <- as.vector(unlist(ww_box_stretch.l))


# --- zebrafish zygotic centred counts --- #
# read domTSS object centred on prim6
  domTSS_egg_p6_shift_fish_prim6_centred.grl <- readRDS("zebrafish/domTSS_egg_p6_shift_fish_prim6_centred_grl.RDS")

# select 10 bp from -34 to -23
  up <- 0
  down <- 12
  range <- c(-up, down)
  win <- up + down

  domTSS_win.gr <- promoters(domTSS_egg_p6_shift_fish_prim6_centred.grl[[1]], up = up, down = down)
  domTSS_win.gr <- shift(domTSS_win.gr, shift = -34*as.integer(ifelse(strand(domTSS_win.gr) =="-", -1, 1)))

# extract sequence 
  domTSS_seq <- getSeq(BSgenome.Drerio.UCSC.danRer7, domTSS_win.gr)
  domTSS.l <- as.list(as.character(domTSS_seq))
  domTSS_char.l <- lapply(domTSS.l, function(x) {
    s2c(x)
  })

# count WW stretches
  word_size <- 12
  ww_box_stretch.l <- lapply(domTSS_char.l, function(x) {
    ww_box_length(x, word_size)
  })

  ww_box_stretch_fish_zygotic <- as.vector(unlist(ww_box_stretch.l))


### --- mouse  somatic (mESC) centred counts--- ###
  domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/domTSS_PGC_oocyte_embryo_vs_mESC_shift_grl.RDS")
  # select 10 bp from -34 to -23
  up <- 0
  down <- 12
  range <- c(-up, down)
  win <- up + down

  domTSS_win.gr <- promoters(domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl[["PN6"]], up = up, down = down)
  domTSS_win.gr <- shift(domTSS_win.gr, shift = -34*as.integer(ifelse(strand(domTSS_win.gr) =="-", -1, 1)))

# extract sequence 
  domTSS_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, domTSS_win.gr)
  domTSS.l <- as.list(as.character(domTSS_seq))
  domTSS_char.l <- lapply(domTSS.l, function(x) {
    s2c(x)
  })

# count WW stretches
  word_size <- 12
  ww_box_stretch.l <- lapply(domTSS_char.l, function(x) {
    ww_box_length(x, word_size)
  })

  ww_box_stretch_mouse_somatic <- as.vector(unlist(ww_box_stretch.l))

# ----plotting distributions
# convert to ggplot df
  ww_box_fish_domTSS_egg.df <- data.frame("counts" = ww_box_stretch_fish, "sample" = "ww_box_fish_domTSS_egg")
  ww_box__mouse_domTSS_egg.df <- data.frame("counts" = ww_box_stretch_mouse, "sample" = "ww_box__mouse_domTSS_egg")
  ww_box_fish_domTSS_som.df <- data.frame("counts" = ww_box_stretch_fish_zygotic, "sample" = "ww_box_fish_domTSS_som")
  ww_box_mouse_domTSS_som.df <- data.frame("counts" = ww_box_stretch_mouse_somatic, "sample" = "ww_box_mouse_domTSS_som")

# combine into one dataframe
  ww_box_counts.df <- rbind(ww_box_fish_domTSS_egg.df, ww_box__mouse_domTSS_egg.df, ww_box_fish_domTSS_som.df, ww_box_mouse_domTSS_som.df)

# change levels
  ww_box_counts.df$sample <-factor(ww_box_counts.df$sample, levels = c("ww_box__mouse_domTSS_egg", "ww_box_fish_domTSS_egg", "ww_box_mouse_domTSS_som", "ww_box_fish_domTSS_som"))

# ---- select mouse distributions - plots egg vs mESCE14 centred---- #
library(magrittr)
library(dplyr)

ww_box_counts_sel_mouse.df <- as_tibble(ww_box_counts.df) %>% filter(sample  == "ww_box__mouse_domTSS_egg" | sample == "ww_box_mouse_domTSS_som") 

# drop unused levels
  ww_box_counts_sel_mouse.df <- droplevels(ww_box_counts_sel_mouse.df)

# plot distribution of domTSS distances 
  col <-c("#440154FF", "firebrick1")

# calculate medians for each distribution
  mouse_med_egg <- median(ww_box_counts_sel_mouse.df %>% filter(sample =="ww_box__mouse_domTSS_egg") %>% dplyr::select(counts) %>% unlist)
  mouse_med_som <- median(ww_box_counts_sel_mouse.df %>% filter(sample =="ww_box_mouse_domTSS_som") %>% dplyr::select(counts) %>% unlist)

  plotA <- ggplot(ww_box_counts_sel_mouse.df, aes(x = counts,  fill = sample, alpha = 0.85)) +
    geom_histogram(aes(x = ww_box_counts_sel_mouse.df$counts , y = ..count../sum(..count..)*100), col = "black", lwd = 0.5, binwidth = 1, position = "identity", alpha = 0.4) +
    scale_fill_manual("sample", values = c("#440154FF", "firebrick1")) +
    theme_light()  +
    theme(text = element_text(size = 20, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          legend.position = "top") +
    labs(y = "%", x = "WW-box length") + 
    scale_x_continuous(limits = c(-1, 13), breaks = seq(0, 13, by = 2)) +
    scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 5)) +
    geom_vline(xintercept = mouse_med_egg, lty = "dashed", lwd = 1, col = col[1]) +
    geom_vline(xintercept = mouse_med_som, lty = "dashed", lwd = 1, col = col[2])

# ---- select zebrafish distributions - plots egg vs prim6---- #
library(magrittr)
library(dplyr)

  ww_box_counts_sel_fish.df <- as_tibble(ww_box_counts.df) %>% filter(sample  == "ww_box_fish_domTSS_egg" | sample == "ww_box_fish_domTSS_som") 

# drop unused levels
  ww_box_counts_sel_fish.df <- droplevels(ww_box_counts_sel_fish.df)

# plot distribution of domTSS distances 
  col <- c("#440154FF", "firebrick1") # colour firebrick corresponds to col scheme firebrick used for mESC

# calculate medians for each distribution
  fish_med_egg <- median(ww_box_counts_sel_fish.df %>% filter(sample =="ww_box_fish_domTSS_egg") %>% dplyr::select(counts) %>% unlist)
  fish_med_som <- median(ww_box_counts_sel_fish.df %>% filter(sample =="ww_box_fish_domTSS_som") %>% dplyr::select(counts) %>% unlist)

  plotB <- ggplot(ww_box_counts_sel_fish.df, aes(x = counts,  fill = sample, alpha = 0.85)) +
    geom_histogram(aes(x =ww_box_counts_sel_fish.df$counts , y = ..count../sum(..count..)*100), col = "black", lwd = 0.5, binwidth = 1, position = "identity", alpha = 0.4) +
    scale_fill_manual("sample", values = c("#440154FF", "firebrick1")) +
    theme_light()  +
    theme(text = element_text(size = 20, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          legend.position = "top") +
    labs(y = "%", x = "WW-box length") + 
    scale_x_continuous(limits = c(-1, 13), breaks = seq(0, 13, by = 2)) +
    scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 5)) +
    geom_vline(xintercept = fish_med_egg, lty = "dashed", lwd = 1, col = col[1]) +
    geom_vline(xintercept = fish_med_som, lty = "dashed", lwd = 1, col = col[2])

# connect to 1 figure 2 panels
  library(cowplot)

pdf("Figures/WW_box_stretch_length_distrib.pdf", width = 10, height = 4)
  plot_grid(plotA, plotB, align = "h", ncol = 2)
dev.off() 
