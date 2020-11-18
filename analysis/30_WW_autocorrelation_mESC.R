#/*==========================================================================#*/
#' ## Extended Data Figure 12M
#/*==========================================================================#*/
# Autocorrelation analyses of WW dinucleotides in mESC broad promoters - all and subsampled

### WW signal for canonical broad promoters - mESC 
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(magrittr)

# load domTSS object
  domTSS_PGCs_oocyte_embryo_anno.grl <- readRDS(file ="../all_reps_PN6_2cell/intermediate_data/domTSS_PGCs_oocyte_embryo_anno_grl.RDS")
  samples <- names(domTSS_PGCs_oocyte_embryo_anno.grl)

# repeating WW extraction on [+50 to +200bp] and no smoothing, binsize = 1
# create windows centered on domTSS 
  up <- 0
  down <- 200
  range <- c(-up, down)
  win <- up + down

# select only mESC
  domTSS_E14_win.grl <- sapply(domTSS_PGCs_oocyte_embryo_anno.grl, function(x) promoters(x, up = up, down = down))

# remove out of bound ranges
  domTSS_E14_win.grl <- sapply(domTSS_E14_win.grl, function(x) x[width(trim(x)) == win])

# extract sequence 
  domTSS_E14_seq.l <- sapply(domTSS_E14_win.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))
  domTSS_E14_seq.l <- list(domTSS_E14_seq.l[[1]])

# use malcolms heatmaps to extract data for metaplot - see dataMetaplot function at end of code
  library(heatmaps)

# plot dinucleotide plots centered on dominant TSS
  pattern <- "WW"
  broad_idx <- domTSS_E14_win.grl[[1]]$interquantile_width >= 9
  sharp_idx <- domTSS_E14_win.grl[[1]]$interquantile_width < 9

# plot metaplots 
# produce heatmaps prior to plotting (no smoothing)
  hm_broad.l <- list()
  hm_sharp.l <- list()

  for(i in 1:length(domTSS_E14_seq.l)) {
    hm_broad.l[[i]] <- PatternHeatmap(domTSS_E14_seq.l[[i]][broad_idx], pattern = pattern, coords = range, label = paste(pattern, ", ", samples[i], sep = ""))
    hm_sharp.l[[i]] <- PatternHeatmap(domTSS_E14_seq.l[[i]][sharp_idx], pattern = pattern, coords = range, label = paste(pattern, ", ", samples[i], sep = ""))
  }

# convert heatmaps to metaplot data - returns df 
  meta_broad.l <- lapply(hm_broad.l, function(x) dataMetaplot(x, binsize = 1))
  meta_sharp.l <- lapply(hm_sharp.l, function(x) dataMetaplot(x, binsize = 1))

# create a list of dataframes (with sharp and broad column)
  all_df.l <- list()

  for (i in 1:length(domTSS_E14_seq.l)) {
    all_df.l[[i]] <- data.frame("x_coord" = meta_broad.l[[i]]$x_coord, 
                                "broad_occurrence" = meta_broad.l[[i]]$occurrence,
                                "sharp_occurrence" = meta_sharp.l[[i]]$occurrence)
  }

# collapse the list into a dataframe
  all_df <- do.call(rbind, all_df.l)
  rownames(all_df) <- 1:nrow(all_df)

# filter to take coordinate from 50 to 200 bp
  all_df <- all_df %>% dplyr::filter(x_coord > 49)

# plot metaplots using ggplot2
  library(ggplot2)
  p <- ggplot(all_df) +
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
    #scale_y_continuous(name = "WW relative frequency", limits = c(0.12, 0.15), seq(0.12, 0.15, by = 0.05)) +  
    scale_x_continuous(name = "Distance to dominant TSS (bp)", limits = c(50, 201), seq(50, 201, by = 10))  

# calculate autocorrelation for each sample
  dsharp <- data.frame("autocorr" = acf(all_df %>% dplyr::pull(sharp_occurrence), plot = F)$acf,
                       "class" = "sharp",
                       "pos" = c(1:22))

  dbroad <- data.frame("autocorr" = acf(all_df %>% dplyr::pull(broad_occurrence), plot = F)$acf,
                       "class" = "broad",
                       "pos" = c(1:22))  

  autocorr_samples_df <- rbind(dsharp, dbroad)

# change plotting levels random vs domTSS
  autocorr_samples_df$class <- factor(autocorr_samples_df$class, levels = c("broad", "sharp"))

## -- autocorrelation plotting --##
# all sequences
  library(ggplot2)
  pA <- ggplot(data = autocorr_samples_df, aes(x = pos,
                                               y = autocorr,
                                               fill = class)) + 
    scale_fill_manual(values = c("royalblue4", "goldenrod2")) +
    geom_bar(stat = "identity", position = position_dodge(), colour = "black", size = 0.25, alpha = 0.6) + 
    geom_line(aes(colour =  class)) +
    scale_colour_manual(values = c("royalblue4", "goldenrod2"), name = NULL) +
    theme_light() +
    theme(text = element_text(size = 20, colour = "black"), 
          axis.text.x = element_text(size = 20, colour = "black"),
          axis.text.y = element_text(size = 20, colour = "black")) 



# select subset of sequences for metaplot - corresponding to numbers in E16.5_F - 560
#broad
  broad.seq.l <- domTSS_E14_seq.l[[1]][broad_idx]

# sample 560 sequences
  set.seed(7)
  broad.seq.random.l <- broad.seq.l[sample(length(broad.seq.l), 560)]

#sharp
  sharp.seq.l <- domTSS_E14_seq.l[[1]][sharp_idx]
  
# sample 560 sequences
  set.seed(7)
  sharp.seq.random.l <- sharp.seq.l[sample(length(sharp.seq.l), 560)]
  
  for(i in 1:length(domTSS_E14_seq.l)) {
    hm_broad.l[[i]] <- PatternHeatmap(broad.seq.random.l, pattern = pattern, coords = range, label = paste(pattern, ", ", samples[i], sep = ""))
    hm_sharp.l[[i]] <- PatternHeatmap(sharp.seq.random.l, pattern = pattern, coords = range, label = paste(pattern, ", ", samples[i], sep = ""))
  }

# convert heatmaps to metaplot data - returns df 
  meta_broad.l <- lapply(hm_broad.l, function(x) dataMetaplot(x, binsize = 1))
  meta_sharp.l <- lapply(hm_sharp.l, function(x) dataMetaplot(x, binsize = 1))

# create a list of dataframes (with sharp and broad column)
  all_df.l <- list()

  for (i in 1:length(domTSS_E14_seq.l)) {
    all_df.l[[i]] <- data.frame("x_coord" = meta_broad.l[[i]]$x_coord, 
                                "broad_occurrence" = meta_broad.l[[i]]$occurrence,
                                "sharp_occurrence" = meta_sharp.l[[i]]$occurrence)
  }


# collapse the list into a dataframe
  all_df <- do.call(rbind, all_df.l)
  rownames(all_df) <- 1:nrow(all_df)

# filter to take coordinate from 50 to 200 bp
  all_df <- all_df %>% dplyr::filter(x_coord > 49)

# calculate autocorrelation for each sample
  dsharp <- data.frame("autocorr" = acf(all_df %>% dplyr::pull(sharp_occurrence), plot = F)$acf,
                       "class" = "sharp",
                       "pos" = c(1:22))

  dbroad <- data.frame("autocorr" = acf(all_df %>% dplyr::pull(broad_occurrence), plot = F)$acf,
                       "class" = "broad",
                       "pos" = c(1:22))  

  autocorr_samples_df <- rbind(dsharp, dbroad)

# change plotting levels random vs domTSS
  autocorr_samples_df$class <- factor(autocorr_samples_df$class, levels = c("broad", "sharp"))

## -- autocorrelation plotting --##
# all sequences
  pB <- ggplot(data = autocorr_samples_df, aes(x = pos,
                                               y = autocorr,
                                               fill = class)) + 
    scale_fill_manual(values = c("royalblue4", "goldenrod2")) +
    geom_bar(stat = "identity", position = position_dodge(), colour = "black", size = 0.25, alpha = 0.6) + 
    geom_line(aes(colour =  class)) +
    scale_colour_manual(values = c("royalblue4", "goldenrod2"), name = NULL) +
    theme_light() +
    theme(text = element_text(size = 20, colour = "black"), 
          axis.text.x = element_text(size = 20, colour = "black"),
          axis.text.y = element_text(size = 20, colour = "black")) 


# connect panels using cowplot
  library(cowplot)
  pdf("Figures/WW_50_200_autocorr_mESC_domTS_sharp_broad.pdf", width = 12, height = 4)
  plot_grid(pA, pB, align = "h", ncol = 2)
  dev.off() 


## ---------- Malcolm's metaplot function from Heatmaps package with plotting removed ---------- ##
#' Plot a Meta-region plot from heatmaps
#'
#' @param hm_list A list of heatmaps
#' @param binsize Integer, size of bins to use in plot
#' @param colors Color to use for each heatmap
#' @param addReferenceLine Logical, add reference line at zero or not
#'
#' This function creates a meta-region plot from 1 or more heatmaps with the same
#' coordinates. A meta-region plot graphs the sum of the signal at each position in
#' each heatmap rather than visualising the signal in two dimensions. Often binning
#' is required to smooth noisy signal.
#'
#' @return invisible(0)
#'
#' @export
#' @importFrom graphics plot mtext legend axis lines
#' @examples
#' data(HeatmapExamples)
#' plotHeatmapMeta(hm, color="steelblue")

dataMetaplot = function(hm_list, binsize=1, colors=gg_col(length(hm_list)), addReferenceLine=FALSE) {
  if (class(hm_list) == "Heatmap") hm_list = list(hm_list) # allow single heatmap argument
  
  if (!length(unique(lapply(hm_list, function(x) x@coords))) == 1)
    stop("heatmaps must have the same coordinates")
  
  
  if (!length(unique(lapply(hm_list, function(x) xm(x)))) == 1)
    stop("heatmaps must have the same xm values")
  
  n_seq = unique(sapply(hm_list, function(x) x@nseq))
  
  if (!length(n_seq) == 1)
    stop("heatmaps must have the same number of sequences")
  
  coords = hm_list[[1]]@coords
  if (binsize != 1) {
    if (!all(xm(hm_list[[1]]) == 1:width(hm_list[[1]]))) {
      stop("cannot set binsize for heatmaps which are already binned/smoothed")
    }
    breaks = seq(0, width(hm_list[[1]]), by=binsize)
    bin_sums = lapply(hm_list, bin_heatmap, breaks=breaks)
    x_coord = breaks[1:(length(breaks)-1)] + binsize/2 + coords[1]
  } else {
    breaks = xm(hm_list[[1]])
    x_coord = breaks + binsize/2 + coords[1]
    bin_sums = lapply(hm_list, function(x) colSums(image(x)))
  }
  scale_factor = n_seq*(width(hm_list[[1]])/length(x_coord))
  occurrence = lapply(bin_sums, function(x) x/scale_factor)
  max_value = max(vapply(occurrence, max, numeric(1)))
  output = data.frame("x_coord" = x_coord, "occurrence" = occurrence[[1]], "bin_sums" = bin_sums[[1]])
  return(output)
}

#' @importFrom stats aggregate
bin_heatmap = function(hm, breaks) {
  partition = data.frame(pos=xm(hm), value=colSums(image(hm)), bin=cut(xm(hm), breaks))
  aggregate(partition$value, sum, by=list(bin=partition$bin))$x
}

#' @importFrom grDevices hcl
gg_col = function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

