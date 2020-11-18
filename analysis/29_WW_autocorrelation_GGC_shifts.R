#/*==========================================================================#*/
#' ## Extended Data Figure 12L, M
#/*==========================================================================#*/
# Autocorrelation analyses of WW dinucleotides 

### Autocorrelation analyses for all WW dinucleotides ---- ###
# repeating WW extraction on [+50 to +200bp] and no smoothing, binsize = 1
# create windows centered on domTSS 
  up <- 0
  down <- 200
  range <- c(-up, down)
  win <- up + down

  domTSS_E14_win.grl <- sapply(domTSS_overlap.l, function(x) promoters(x, up = up, down = down))

# remove out of bound ranges
  domTSS_E14_win.grl <- sapply(domTSS_E14_win.grl, function(x) x[width(trim(x)) == win])

# extract sequence 
  domTSS_E14_seq.l <- sapply(domTSS_E14_win.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))

# use malcolms heatmaps sum data it into a metaplot..- plot using ggplot2
  library(heatmaps)

# plot dinucleotide plots centered on dominant TSS
  pattern <- "WW"
  hm.l <- list()
  hm_smoothed.l <- list()

  names(sort.broad.l) <- names(domTSS_E14_win.grl)

  # sort index - by distance
  sort.l <- lapply(domTSS_E14_win.grl, function(x) {
    return(order(x$domTSS_dist, decreasing = T))
  })

# plot metaplots 
  # produce heatmaps prior to plotting (no smoothing) - see function at end of code
    hm_broad.l <- list()

    for(i in 1:length(domTSS_E14_seq.l)) {
      hm_broad.l[[i]] <- PatternHeatmap(domTSS_E14_seq.l[[i]], pattern = pattern, coords = range, label = paste(pattern, ", ", samples[i], sep = ""))
    }
    names(hm_broad.l) <- names(domTSS_E14_seq.l)

  # convert heatmaps to metaplot data - returns df 
    meta_broad.l <- lapply(hm_broad.l, function(x) dataMetaplot(x, binsize = 1))

# create a list of dataframes (with sharp and broad column)
  all_df.l <- list()

  samples <- names(domTSS_E14_seq.l)

  for (i in 1:length(samples)) {
    all_df.l[[i]] <- data.frame("x_coord" = meta_broad.l[[i]]$x_coord, 
                                "broad_occurrence" = meta_broad.l[[i]]$occurrence,
                                "sample" = rep(samples[i], times = length(meta_broad.l[[1]]$occurrence)))
  }
  names(all_df.l) <- samples


# remove tbp2 KO
  all_df.l <- all_df.l[-c(14, 16)]

# collapse the list into a dataframe
  all_df <- do.call(rbind, all_df.l)
  rownames(all_df) <- 1:nrow(all_df)

# filter to take coordinate from 50 to 200 bp
  all_df <- all_df %>% filter(x_coord > 49)


# -- random TSS windows -- #
# # create windows centered on random TSS from that tag cluster +50 to +200
  up <- 0
  down <- 200
  range <- c(-up, down)
  win <- up + down

  randomTSS_E14_win.grl <- sapply(random.ctss.coord.l, function(x) promoters(x, up = up, down = down))

# remove out of bound ranges
  randomTSS_E14_win.grl <- sapply(randomTSS_E14_win.grl, function(x) x[width(trim(x)) == win])

# check the length of the each range in the list
  sapply(randomTSS_E14_win.grl, length)

# extract sequence 
  randomTSS_E14_seq.l <- sapply(randomTSS_E14_win.grl, function(x) getSeq(BSgenome.Mmusculus.UCSC.mm10, x))

# exclude TBP2KO
  randomTSS_E14_seq.l <- randomTSS_E14_seq.l[-c(14, 16)]

# use malcolms heatmaps without plotting
  library(heatmaps)

# plot dinucleotide plots centered on dominant TSS
  pattern <- "CC"
  hm.l <- list()
  hm_smoothed.l <- list()

# plot metaplots 
# produce heatmaps prior to plotting (no smoothing)
  hm_broad.l <- list()
  
  for(i in 1:length(randomTSS_E14_seq.l)) {
    hm_broad.l[[i]] <- PatternHeatmap(randomTSS_E14_seq.l[[i]], pattern = pattern, coords = range, label = paste(pattern, ", ", samples[i], sep = ""))
  }
  names(hm_broad.l) <- names(randomTSS_E14_seq.l)

# convert heatmaps to metaplot data - returns df 
  meta_broad.l <- lapply(hm_broad.l, function(x) dataMetaplot(x, binsize = 1))

# create a list of dataframes (with sharp and broad column)
  random_df.all <- list()

  samples <- names(randomTSS_E14_seq.l)

  for (i in 1:length(samples)) {
    random_df.all[[i]] <- data.frame("x_coord" = meta_broad.l[[i]]$x_coord, 
                                     "broad_occurrence" = meta_broad.l[[i]]$occurrence,
                                     "sample" = rep(samples[i], times = length(meta_broad.l[[1]]$occurrence)),
                                     "class" = rep("random", times = length(meta_broad.l[[1]]$occurrence)))
  }
  names(random_df.all) <- samples

# collapse the list into a dataframe
  random_df <- do.call(rbind, random_df.all)
  rownames(random_df) <- 1:nrow(random_df)

# filter random df to include coordinates 50 bp to 200 bp
  random_df <- random_df %>% filter(x_coord > 49)

  all_df$class = "domTSS"

# combine random and normal TSS data
  all_combined.gg <- rbind(all_df, random_df)


# calculate autocorrelation for each sample
  samples <- levels(all_combined.gg$sample)
  autocorr_samples.l <- lapply(samples, function(x) {
    dr <- data.frame("autocorr" = acf(all_combined.gg %>% filter(sample == x, class == "random") %>% pull(broad_occurrence), plot = F)$acf,
                     "sample" = x,
                     "class" = "random",
                     "pos" = c(1:22))
    dd <- data.frame("autocorr" = acf(all_combined.gg %>% filter(sample == x, class == "domTSS") %>% pull(broad_occurrence), plot = F)$acf,
                     "sample" = x,
                     "class" = "domTSS",
                     "pos" = c(1:22))  
    return(rbind(dr,dd))
  })

# collapse the autocorrelationlist into a dataframe
  autocorr_samples_df <- do.call(rbind, autocorr_samples.l)

# change plotting levels random vs domTSS
  autocorr_samples_df$class <- factor(autocorr_samples_df$class, levels = c("domTSS", "random"))

## -- autocorrelation plotting --##
  p <- ggplot(data = autocorr_samples_df, aes(x = pos,
                                              y = autocorr,
                                              fill = class)) + 
    scale_fill_manual(values = c("royalblue4", "goldenrod2"), name = NULL) +
    geom_bar(stat = "identity", position = position_dodge(), colour = "black", size = 0.25, alpha = 0.6) + 
    geom_line(aes(colour =  class)) +
    scale_colour_manual(values = c("royalblue4", "goldenrod2"), name = NULL) +
    theme_light() +
    theme(text = element_text(size = 20, colour = "black"), 
          axis.text.x = element_text(size = 20, colour = "black"),
          axis.text.y = element_text(size = 20, colour = "black")) +
    guides(fill = FALSE)
  
  p <- p + facet_wrap(~ sample, ncol = 4)
  
  pdf("Figures/WW_50_200_autocorr_all_domTSS_random.pdf", height = 12, width = 12)
    print(p)
  dev.off()


# print selected samples - E16_5F and PN6
  selection.gg <- autocorr_samples_df %>% filter(sample == "E16_5F" | sample == "PN6")
  selection.gg$sample <- droplevels(selection.gg$sample)
  
  p <- ggplot(data = selection.gg, aes(x = pos,
                                       y = autocorr,
                                       fill = class)) + 
    scale_fill_manual(values = c("royalblue4", "goldenrod2"), name = NULL) +
    geom_bar(stat = "identity", position = position_dodge(), colour = "black", size = 0.25, alpha = 0.6) + 
    geom_line(aes(colour =  class)) +
    scale_colour_manual(values = c("royalblue4", "goldenrod2"), name = NULL) +
    theme_light() +
    theme(text = element_text(size = 20, colour = "black"), 
          axis.text.x = element_text(size = 20, colour = "black"),
          axis.text.y = element_text(size = 20, colour = "black")) +
    guides(fill = FALSE)
  
  p <- p + facet_wrap(~ sample, ncol = 2)

  pdf("Figures/WW_50_200_autocorr_select_domTSS_random.pdf", height = 4, width = 11)
    print(p)
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

