#/*==========================================================================#*/
#' ## Extended Data Figure 15D,E
#/*==========================================================================#*/
# WW autocorrelation in dominant CTSS SOM classes

# repeating WW extraction on [+50 to +200bp] and no smoothing, binsize = 1
# create windows centered on domTSS 
  up <- 0
  down <- 200
  range <- c(-up, down)
  win <- up + down
  
  domTSS_win.grl <-lapply(CTSSs_SOMclass_sample_filtr.grl, function(grl) {
    
    soms_win.grl <- lapply(grl, function(gr) {
      domTSS_win.gr <- promoters(gr, up = up, down = down)
      
      # remove out of bound ranges
      domTSS_win.gr <- domTSS_win.gr[width(trim(domTSS_win.gr)) == win]
      return(domTSS_win.gr)
      
      # drop mitochondrial sequences
      domTSS_win.gr <- dropSeqlevels(domTSS_win.gr, value = "chrM", pruning.mode = "coarse")
    })
    return(soms_win.grl)
  }) 


# extract sequence 
  domTSS_seq.l <- lapply(domTSS_win.grl, function(grl) {
    soms_seq.l <- lapply(grl, function(gr) {
      seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, gr)
      return(seq)
    })
    return(soms_seq.l)
  })

# use malcolms heatmaps to create smoothed pattern heatmaps and sum it into a metaplot..- plot using ggplot2
  library(heatmaps)

# plot dinucleotide plots centered on dominant TSS
  pattern <- "WW"
  hm.l <- list()
  hm_smoothed.l <- list()

# plot metaplots 
# produce heatmaps prior to plotting (no smoothing)
  hm_broad.l <- list()
  
  for(i in 1:length(domTSS_seq.l)) {
    for(j in 1:length(domTSS_seq.l[[i]])) {
      hm_broad.l[[j]] <- PatternHeatmap(domTSS_seq.l[[i]][[j]], pattern = pattern, coords = range, label = paste(pattern, ", ", names(domTSS_seq.l)[j], names(domTSS_seq.l[[i]])[j], "_", sep = ""))
      
    }
    names(hm_broad.l) <- names(domTSS_seq.l[[i]])
  
# convert heatmaps to metaplot data - returns df - see function at the end of script
  meta_broad.l <- lapply(hm_broad.l, function(x) metaplot_data(x, up = 0, down = 200))
  
# create a list of dataframes (with sharp and broad column)
  all_df.l <- list()
  
  samples <- names(domTSS_seq.l[[i]])
  
  for (j in 1:length(samples)) {
    all_df.l[[j]] <- data.frame("x_coord" = meta_broad.l[[j]]$coord, 
                                "broad_occurrence" = meta_broad.l[[j]]$df_smooth.l,
                                "sample" = rep(samples[j], times = length(meta_broad.l[[1]]$df_smooth.l)))
  }
  names(all_df.l) <- samples
  
  
  # collapse the list into a dataframe
    all_df <- do.call(rbind, all_df.l)
    rownames(all_df) <- 1:nrow(all_df)
  
  
  # filter to take coordinate from 50 to 200 bp
    all_df <- all_df %>% filter(x_coord > 49)
  
  # calculate autocorrelation for each sample
    samples <- levels(all_df$sample)
    autocorr_samples.l <- lapply(samples, function(x) {
      
      dd <- data.frame("autocorr" = acf(all_df %>% filter(sample == x) %>% pull(broad_occurrence), plot = F)$acf,
                       "sample" = x,
                       "class" = "domTSS",
                       "pos" = c(1:22))  
      return(rbind(dd))
    })
  
  # collapse the autocorrelationlist into a dataframe
    autocorr_samples_df <- do.call(rbind, autocorr_samples.l)
  
  ## -- autocorrelation plotting --##
    p <- ggplot(data = autocorr_samples_df, aes(x = pos,
                                                y = autocorr)) +
      geom_bar(stat = "identity", position = position_dodge(), colour = "black", size = 0.25, alpha = 0.4, fill = "royalblue4") + 
      geom_line(inherit.aes = T, colour =  "royalblue4") +
      theme_light() +
      theme(text = element_text(size = 20, colour = "black"), 
            axis.text.x = element_text(size = 20, colour = "black"),
            axis.text.y = element_text(size = 20, colour = "black")) +
      guides(fill = FALSE)
    
    p <- p + facet_wrap(~ sample, ncol = 5)
    
    pdf(paste0("Figures/seqlogo_domCTSS_som/WW_50_200_autocorr_domTSS_", names(domTSS_seq.l)[i], ".pdf"), height = 12, width = 12)
    print(p)
    dev.off()
    
  }

  # ---- metaplot_data function ---- #
  metaplot_data <- function(hm.l, motif = FALSE, up = 500,  down = 500, smooth_span = 3, smooth_step = 1) {
    
    # path - sets where to save the image and image name
    # hm.l is a heatmap object (not a list) - Malcolm's heatmap object with extractable image
    # classes is a factor  class that divides a metaplot
    
    # subset matrix in classes
    matrix_l <- image(hm.l)
    
    # metaplot calculation and smoothing
    library(zoo)
    matrix_sum.l <-  colSums(matrix_l)/nrow(matrix_l)
    df_smooth.l <- rollapply(data = matrix_sum.l, FUN = mean, width = smooth_span, by = smooth_step)
    
    df_smooth.df <- as.data.frame(df_smooth.l)
    
    # add coordinate column
    if (motif == FALSE) {
      df_smooth.df$coord <- rollapply(data = -up:(down - 1), FUN = mean, width = smooth_span, by = smooth_step)
    } else {
      df_smooth.df$coord <- rollapply(data = -up:(down - (up + down - nrow(df_smooth.df) - 1)), FUN = mean, width = smooth_span, by = smooth_step)
    }
    
    # format dataframe to ggplot acceptable
    library(tidyr)
    return(df_smooth.df)
  }