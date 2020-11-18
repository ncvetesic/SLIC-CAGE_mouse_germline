#/*==========================================================================#*/
#' ## Extended Data Figure 15B,C
#/*==========================================================================#*/
# WW periodicity in dominant CTSS SOM classes

library(ggplot2)

# load SOM divided dominant CTSS object
  domCTSSs_SOM25_PGC_F_oocyte_promOnly.grl <- readRDS("intermediate_data/domCTSSs_SOM25_PGC_F_oocyte_promOnly.grl.RDS")

# add SOM class column to each GRanges list element  
  CTSSs_SOMclass.grl <- GRangesList(domCTSSs_SOM25_PGC_F_oocyte_promOnly.grl)
  CTSSs_SOMclass.grl <- lapply(1:length(CTSSs_SOMclass.grl), function(x) {
    CTSSs_SOMclass.grl[[x]]$SOM_class <- names(CTSSs_SOMclass.grl)[x]
    return(CTSSs_SOMclass.grl[[x]])
  })

# unlist or concatenate to one GRanges
  CTSSs_SOMclass.gr <- do.call(c, unlist(CTSSs_SOMclass.grl))

# separate into grangeslist by samples
  samples <- colnames(mcols(CTSSs_SOMclass.gr))[-c(1:3, 11, 12)]

  CTSSs_SOMclass_samples.grl <- lapply(1:length(samples), function(x) {
    gr <- CTSSs_SOMclass.gr
    mcols(gr) <- NULL
    gr$tpm <- mcols(CTSSs_SOMclass.gr)[, samples[x]]
    gr$SOM_class <- mcols(CTSSs_SOMclass.gr)$SOM_class
    gr$CTSSid  <- mcols(CTSSs_SOMclass.gr)$CTSSid 
    return(gr)
  })

# add sample names
  names(CTSSs_SOMclass_samples.grl) <- samples

# filter to include only CTSSs with at least 1 tpm - using 0.5 TPM filter gives about 110,000 vs 80,000 CTSSs
  CTSSs_SOMclass_samples_filtr.grl <- lapply(CTSSs_SOMclass_samples.grl, function(x) {
    gr <- x[x$tpm >= 0.5]
    return(gr)
  })

# add sample column and unlist
  CTSSs_SOMclass_samples_filtr.grl <- lapply(1:length(CTSSs_SOMclass_samples_filtr.grl), function(x){
    gr <- CTSSs_SOMclass_samples_filtr.grl[[x]]
    gr$sample <- names(CTSSs_SOMclass_samples_filtr.grl)[x]
    return(gr)
  })


# unlist or concatenate to one GRanges
  CTSSs_SOMclass_sample_filtr.gr <- do.call(c, unlist(CTSSs_SOMclass_samples_filtr.grl))

# divide into per sample and then per som.
  samples <- unique(CTSSs_SOMclass_sample_filtr.gr$sample)
  soms <- unique(CTSSs_SOMclass_sample_filtr.gr$SOM_class)

  CTSSs_SOMclass_sample_filtr.grl <- lapply(1:length(samples), function(x) {
    gr <- CTSSs_SOMclass_sample_filtr.gr[CTSSs_SOMclass_sample_filtr.gr$sample == samples[x]]
    gr.l <- list()
    gr.l <- lapply(1:length(soms), function(x){
      gr2 <- gr[gr$SOM_class == soms[x]]
    })
    names(gr.l) <- soms
    return(gr.l)
  })

# name lists with sample names
  names(CTSSs_SOMclass_sample_filtr.grl) <- samples

# need 2 loops - one for sampples 1 for soms - should plot each seqlogo separately as i wont use all
#--- selecting W-box regions - centered on maternal TSS ---#
# create windows for plotting
  up <- 100
  down <- 300
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
  
# convert heatmaps to metaplot data - returns df - see function at the end of code
  meta_broad.l <- lapply(hm_broad.l, function(x) metaplot_data(x, up = 100, down = 300))
  
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
  
  # collapse df to ggplot compatible df
    library(tidyr)
  
  # plot metaplots using ggplot2
    col <- rep("black", times = 25)
    names(col) <- samples
  
    library(ggplot2)
    p <- ggplot(all_df) +
      geom_line(lty = 1, aes(x = x_coord, y = broad_occurrence, colour = sample), size = 0.6) +
      scale_color_manual(values = col, 
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
      scale_y_continuous(name = "Relative frequency", limits = c(0, 0.5), breaks = c(0, 0.25, 0.5)) +
      scale_x_continuous(name = "Distance to dominant TSS (bp)", limits = c(-100, 300), seq(-100, 300, by = 100))
    
    p <- p + facet_wrap(~ sample, ncol = 5)
    
    pdf(paste0("Figures/seqlogo_domCTSS_som/WW_periodicity_", names(domTSS_seq.l)[i], ".pdf"), height = 12, width = 16)
      print(p)
    dev.off()
    
  
# zoomed metaplots for the inset = only broad pattern occurence
  library(ggplot2)
  p <- ggplot(all_df, aes(x = x_coord, y = broad_occurrence, colour = sample)) +
    geom_line(lty = 1, size = 0.6, legend = FALSE) +
    scale_color_manual(values = col, 
                       name = NULL)  +
    theme_light()   +
    theme(text = element_text(size = 20, colour = "black"),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          legend.position = "right",
          legend.title = element_blank()) +
    guides(colour = guide_legend(reverse=F)) +
    coord_cartesian(xlim = c(50, 200), ylim = c(0.05, 0.2))
  p <- p + facet_wrap(~ sample, ncol = 5)
  
  pdf(paste0("Figures/seqlogo_domCTSS_som/WW_periodicity_zoom_", names(domTSS_seq.l)[i], ".pdf"), height = 12, width = 16)   
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
```