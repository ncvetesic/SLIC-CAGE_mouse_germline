#/*==========================================================================#*/
#' ## Extended Data Figure 12I
#/*==========================================================================#*/
# Correlation of tag cluster signal spread width in E16.5F vs mESC shifting promoters and its CpG content

# bigwig files from Nome-seq paper are prepared using Rtracklayer + liftover from mm9 to mm10 - see example below
library(rtracklayer)
  PGCE16_5F_nome_seq_WCG_rep1.gr <- import("/Volumes/Lucille_data/Mouse_PGC_oocyte_embryo_data/nome_seq/GSM2098127_NOMe-seq_Mouse_E16.5_female_PGC_rep1.WCG.bw")
  PGCE16_5F_nome_seq_WCG_rep2.gr <- import("/Volumes/Lucille_data/Mouse_PGC_oocyte_embryo_data/nome_seq/GSM2098128_NOMe-seq_Mouse_E16.5_female_PGC_rep2.WCG.bw")
  PGCE16_5F_nome_seq_WCG_rep3.gr <- import("/Volumes/Lucille_data/Mouse_PGC_oocyte_embryo_data/nome_seq/GSM2098129_NOMe-seq_Mouse_E16.5_female_PGC_rep3.WCG.bw")

# liftover to mm10
  mm9ToMm10.over.chain<- import.chain("mm9ToMm10.over.chain")
  PGCE16_5F_nome_seq_WCG_rep1_mm10.gr <- unlist(liftOver(PGCE16_5F_nome_seq_WCG_rep1.gr, chain = mm9ToMm10.over.chain))
  PGCE16_5F_nome_seq_WCG_rep2_mm10.gr <- unlist(liftOver(PGCE16_5F_nome_seq_WCG_rep2.gr, chain = mm9ToMm10.over.chain))
  PGCE16_5F_nome_seq_WCG_rep3_mm10.gr <- unlist(liftOver(PGCE16_5F_nome_seq_WCG_rep3.gr, chain = mm9ToMm10.over.chain))

  seqinfo(PGCE16_5F_nome_seq_WCG_rep1_mm10.gr) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(PGCE16_5F_nome_seq_WCG_rep1_mm10.gr)]
  seqinfo(PGCE16_5F_nome_seq_WCG_rep2_mm10.gr) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(PGCE16_5F_nome_seq_WCG_rep2_mm10.gr)]
  seqinfo(PGCE16_5F_nome_seq_WCG_rep3_mm10.gr) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(PGCE16_5F_nome_seq_WCG_rep3_mm10.gr)]

### Spread in E16_5F vs Nome-seq
# select regions of spreading - both ranges have to be converted into a GRangesList as each range difference can produce more regions of difference - hence it has to be a list for each pair we are comparing
# broadened E16_5F
  mESC_vs_E16_5F_broader.gr <- mESC_vs_E16_5F.gr[mESC_vs_E16_5F.gr$interquantile_width < mESC_vs_E16_5F.gr$interquantile_width_E16_5F]

# narrower E16_5F
  mESC_vs_E16_5F_narrower.gr <- mESC_vs_E16_5F.gr[mESC_vs_E16_5F.gr$interquantile_width > mESC_vs_E16_5F.gr$interquantile_width_E16_5F]

# define broad and narrow E16_5 shifting as mESC or E16.6 tcs - start and end ranges - mESC are already defined
# broader E16_5F shifting promoters
  E16_5F_TCs_mESC_vs_E16_5F_broader.gr <- GRanges(seqnames = seqnames(mESC_vs_E16_5F_broader.gr),
                                                  ranges = IRanges(start = mESC_vs_E16_5F_broader.gr$start_E16_5F, end = mESC_vs_E16_5F_broader.gr$end_E16_5F),
                                                  strand = strand(mESC_vs_E16_5F_broader.gr))
  mcols(E16_5F_TCs_mESC_vs_E16_5F_broader.gr) <- mcols(mESC_vs_E16_5F_broader.gr)
  seqinfo(E16_5F_TCs_mESC_vs_E16_5F_broader.gr) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(E16_5F_TCs_mESC_vs_E16_5F_broader.gr)]


# narrower E16_5F shifting promoters
  E16_5F_TCs_mESC_vs_E16_5F_narrower.gr <- GRanges(seqnames = seqnames(mESC_vs_E16_5F_narrower.gr),
                                                   ranges = IRanges(start = mESC_vs_E16_5F_narrower.gr$start_E16_5F, end = mESC_vs_E16_5F_narrower.gr$end_E16_5F),
                                                   strand = strand(mESC_vs_E16_5F_narrower.gr))
  mcols(E16_5F_TCs_mESC_vs_E16_5F_narrower.gr) <- mcols(mESC_vs_E16_5F_narrower.gr)
  seqinfo(E16_5F_TCs_mESC_vs_E16_5F_narrower.gr) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(E16_5F_TCs_mESC_vs_E16_5F_narrower.gr)]
  
  mESC_vs_E16_5F_broader.gr <- as(mESC_vs_E16_5F_broader.gr, "GRangesList")
  E16_5F_TCs_mESC_vs_E16_5F_broader.gr <- as(E16_5F_TCs_mESC_vs_E16_5F_broader.gr, "GRangesList")
  
  E16_5_vs_mESC_setdiff.grl <- setdiff(E16_5F_TCs_mESC_vs_E16_5F_broader.gr, mESC_vs_E16_5F_broader.gr)

# unlist
  E16_5_vs_mESC_setdiff.gr <- unlist(E16_5_vs_mESC_setdiff.grl)

# width
  E16_5_vs_mESC_setdiff_width <- width(E16_5_vs_mESC_setdiff.gr)

# centred on mESC - broader shifting promoters
# --- overlap mESC shifting tag clusters signal with methylation
  overlap_E16_5F_WCG.l <- lapply(PGCE16_5F_nome_seq_WCG_mm10.grl, function(x) findOverlapPairs(E16_5_vs_mESC_setdiff.gr, x))
  overlap_E11_5F_WCG.l <- lapply(PGCE11_5F_nome_seq_WCG_mm10.grl, function(x) findOverlapPairs(E16_5_vs_mESC_setdiff.gr, x))
  overlap_soma_E16_5F_WCG.l <- lapply(somaE16_5F_nome_seq_WCG_mm10.grl, function(x) findOverlapPairs(E16_5_vs_mESC_setdiff.gr, x))
  overlap_E13_5F_WCG.l <- lapply(PGCE13_5F_nome_seq_WCG_mm10.grl, function(x) findOverlapPairs(E16_5_vs_mESC_setdiff.gr, x))
  overlap_E16_5M_WCG.l <- lapply(PGCE16_5M_nome_seq_WCG_mm10.grl, function(x) findOverlapPairs(E16_5_vs_mESC_setdiff.gr, x))
  overlap_E13_5M_WCG.l <- lapply(PGCE13_5M_nome_seq_WCG_mm10.grl, function(x) findOverlapPairs(E16_5_vs_mESC_setdiff.gr, x))

# extract methylation signal
# E16_5F methylation overlapped mESC shifting prom (mESC vs E16.5F) - 1439, 1465, 2401 out of 825 ranges
  overlap_E16_5F_methyl_signal.dfl <- lapply(1:length(overlap_E16_5F_WCG.l), function(x) {
    df <- data.frame(sample = names(overlap_E16_5F_WCG.l)[x], score = overlap_E16_5F_WCG.l[[x]]@second$score) 
    return(df)
  })

# E11_5 methylation overlapped mESC shifting prom (mESC vs E16.5F) - 2442,  2513 out of 825 ranges
  overlap_E11_5F_methyl_signal.dfl <- lapply(1:length(overlap_E11_5F_WCG.l), function(x) {
    df <- data.frame(sample = names(overlap_E11_5F_WCG.l)[x], score = overlap_E11_5F_WCG.l[[x]]@second$score) 
    return(df)
  })

# E16_5F soma - methylation overlapped mESC shifting prom (mESC vs E16.5F) 1459, 1251, 1624  out of 825
  overlap_somaE16_5F_methyl_signal.dfl <- lapply(1:length(overlap_soma_E16_5F_WCG.l), function(x) {
    df <- data.frame(sample = names(overlap_soma_E16_5F_WCG.l)[x], score = overlap_soma_E16_5F_WCG.l[[x]]@second$score) 
    return(df)
  })


# E13_5F soma - methylation overlapped mESC shifting prom (mESC vs E16.5F) 1796, 2435, 2330  out of 825
  overlap_E13_5F_methyl_signal.dfl <- lapply(1:length(overlap_E13_5F_WCG.l), function(x) {
    df <- data.frame(sample = names(overlap_E13_5F_WCG.l)[x], score = overlap_E13_5F_WCG.l[[x]]@second$score) 
    return(df)
  })


# E16_5M  - methylation overlapped mESC shifting prom (mESC vs E16.5F) 1796, 2435, 2330  out of 825
  overlap_E16_5M_methyl_signal.dfl <- lapply(1:length(overlap_E16_5M_WCG.l), function(x) {
    df <- data.frame(sample = names(overlap_E16_5M_WCG.l)[x], score = overlap_E16_5M_WCG.l[[x]]@second$score) 
    return(df)
  })

# E13_5M - methylation overlapped mESC shifting prom (mESC vs E16.5F) 1796, 2435, 2330  out of 825
  overlap_E13_5M_methyl_signal.dfl <- lapply(1:length(overlap_E13_5M_WCG.l), function(x) {
    df <- data.frame(sample = names(overlap_E13_5M_WCG.l)[x], score = overlap_E13_5M_WCG.l[[x]]@second$score) 
    return(df)
  })

# bind into one dataframe
  overlap_E16_5F_methyl_signal.df <- do.call(rbind, overlap_E16_5F_methyl_signal.dfl)
  overlap_E16_5F_methyl_signal.df$sample <- factor(overlap_E16_5F_methyl_signal.df$sample, levels = names(overlap_E16_5F_WCG.l))
  
  overlap_E11_5F_methyl_signal.df <- do.call(rbind, overlap_E11_5F_methyl_signal.dfl)
  overlap_E11_5F_methyl_signal.df$sample <- factor(overlap_E11_5F_methyl_signal.df$sample, levels = names(overlap_E11_5F_WCG.l))
  
  overlap_somaE16_5F_methyl_signal.df <- do.call(rbind, overlap_somaE16_5F_methyl_signal.dfl)
  overlap_somaE16_5F_methyl_signal.df$sample <- factor(overlap_somaE16_5F_methyl_signal.df$sample, levels = names(overlap_soma_E16_5F_WCG.l))
  
  overlap_E13_5F_methyl_signal.df <- do.call(rbind, overlap_E13_5F_methyl_signal.dfl)
  overlap_E13_5F_methyl_signal.df$sample <- factor(overlap_E13_5F_methyl_signal.df$sample, levels = names(overlap_E13_5F_WCG.l))
  
  overlap_E16_5M_methyl_signal.df <- do.call(rbind, overlap_E16_5M_methyl_signal.dfl)
  overlap_E16_5M_methyl_signal.df$sample <- factor(overlap_E16_5M_methyl_signal.df$sample, levels = names(overlap_E16_5M_WCG.l))
  
  overlap_E13_5M_methyl_signal.df <- do.call(rbind, overlap_E13_5M_methyl_signal.dfl)
  overlap_E13_5M_methyl_signal.df$sample <- factor(overlap_E13_5M_methyl_signal.df$sample, levels = names(overlap_E13_5M_WCG.l))
  
# bind into one dataframe
  overlap_methyl_signal_broader.df <- rbind(overlap_E11_5F_methyl_signal.df, overlap_E13_5F_methyl_signal.df, overlap_E13_5M_methyl_signal.df, overlap_E16_5F_methyl_signal.df, overlap_E16_5M_methyl_signal.df, overlap_somaE16_5F_methyl_signal.df)
  

# change levels
  overlap_methyl_signal_broader.df$sample <- factor(overlap_methyl_signal_broader.df$sample, levels = c("PGCE11_5F_WCG_rep1", "PGCE11_5F_WCG_rep2", 
                                                                                                        "PGCE13_5F_WCG_rep1", "PGCE13_5F_WCG_rep2", "PGCE13_5F_WCG_rep3",
                                                                                                        "PGCE13_5M_WCG_rep1", "PGCE13_5M_WCG_rep2", "PGCE13_5M_WCG_rep3",
                                                                                                        "PGCE16_5F_WCG_rep1", "PGCE16_5F_WCG_rep2", "PGCE16_5F_WCG_rep3",
                                                                                                        "PGCE16_5M_WCG_rep1", "PGCE16_5M_WCG_rep2",
                                                                                                        "soma_E16_5F_WCG_rep1", "soma_E16_5F_WCG_rep2", "soma_E16_5F_WCG_rep3"))
## -- plotting methylation signal -- ##
# violin plot
  samples <- c("PGC E11.5 r1", "PGC E11.5 r2", 
               "PGC E13.5F r1", "PGC E13.5F r2", "PGC E13.5F r2", 
               "PGC E13.5M r1", "PGC E13.5M r2", "PGC E13.5M r2", 
               "PGC E16.5F r1", "PGC E16.5F r2", "PGC E16.5F r3", 
               "PGC E16.5M r1", "PGC E16.5M r2",
               "soma E16.5F r1", "soma E16.5F r2", "soma E16.5F r3")

  p <- ggplot(overlap_methyl_signal_broader.df, aes(y = score, x = sample, fill = sample), alpha = 0.6) +
        geom_violin(trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75), size = 0.5) +
        scale_fill_manual(values = c("#D1E11CFF", "#D1E11CFF",
                                     "#8BD646FF", "#8BD646FF", "#8BD646FF",
                                     "#75D054FF", "#75D054FF", "#75D054FF",
                                     "#3FBC73FF", "#3FBC73FF", "#3FBC73FF", 
                                     "#31B57BFF", "#31B57BFF",
                                     "#31B57BFF" , "#31B57BFF" , "#31B57BFF")) +
        scale_x_discrete(labels = samples) +
        theme_light() +
        theme(text = element_text(size = 18),
              axis.text.x = element_text(size = 18, angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 18),
              axis.title.x = element_text(size = 18),
              axis.title.y = element_text(size = 18)) +
        labs( x = "Samples") +
        guides(fill = F) +
        coord_cartesian(ylim = c(-0.2, 0.2))

    pdf("Figures/CpG_methylation_E16_5F_spread_diff_broader.pdf", width = 8, height = 6)
      print(p)
    dev.off()  

# select regions of spreading - both ranges have to be converted into a GRangesList as each range difference can produce more regions of difference - hence it has to be a list for each pair we are comparing
# broadened E16_5F
  mESC_vs_E16_5F_broader.gr <- mESC_vs_E16_5F.gr[mESC_vs_E16_5F.gr$interquantile_width < mESC_vs_E16_5F.gr$interquantile_width_E16_5F]

# narrower E16_5F
  mESC_vs_E16_5F_narrower.gr <- mESC_vs_E16_5F.gr[mESC_vs_E16_5F.gr$interquantile_width > mESC_vs_E16_5F.gr$interquantile_width_E16_5F]


  mESC_vs_E16_5F_narrower.gr <- as(mESC_vs_E16_5F_narrower.gr, "GRangesList")
  E16_5F_TCs_mESC_vs_E16_5F_narrower.gr <- as(E16_5F_TCs_mESC_vs_E16_5F_narrower.gr, "GRangesList")
  
  E16_5_vs_mESC_setdiff.grl <- setdiff(E16_5F_TCs_mESC_vs_E16_5F_narrower.gr, mESC_vs_E16_5F_narrower.gr)

# unlist
  E16_5_vs_mESC_setdiff.gr <- unlist(E16_5_vs_mESC_setdiff.grl)

# width
  E16_5_vs_mESC_setdiff_width <- width(E16_5_vs_mESC_setdiff.gr)

# centred on mESC - narrower shifting promoters
# --- overlap mESC shifting tag clusters signal with methylation
  overlap_E16_5F_WCG.l <- lapply(PGCE16_5F_nome_seq_WCG_mm10.grl, function(x) findOverlapPairs(E16_5_vs_mESC_setdiff.gr, x))
  overlap_E11_5F_WCG.l <- lapply(PGCE11_5F_nome_seq_WCG_mm10.grl, function(x) findOverlapPairs(E16_5_vs_mESC_setdiff.gr, x))
  overlap_soma_E16_5F_WCG.l <- lapply(somaE16_5F_nome_seq_WCG_mm10.grl, function(x) findOverlapPairs(E16_5_vs_mESC_setdiff.gr, x))
  overlap_E13_5F_WCG.l <- lapply(PGCE13_5F_nome_seq_WCG_mm10.grl, function(x) findOverlapPairs(E16_5_vs_mESC_setdiff.gr, x))
  overlap_E16_5M_WCG.l <- lapply(PGCE16_5M_nome_seq_WCG_mm10.grl, function(x) findOverlapPairs(E16_5_vs_mESC_setdiff.gr, x))
  overlap_E13_5M_WCG.l <- lapply(PGCE13_5M_nome_seq_WCG_mm10.grl, function(x) findOverlapPairs(E16_5_vs_mESC_setdiff.gr, x))

# extract methylation signal
# E16_5F methylation overlapped mESC shifting prom (mESC vs E16.5F) - 1439, 1465, 2401 out of 825 ranges
  overlap_E16_5F_methyl_signal.dfl <- lapply(1:length(overlap_E16_5F_WCG.l), function(x) {
    df <- data.frame(sample = names(overlap_E16_5F_WCG.l)[x], score = overlap_E16_5F_WCG.l[[x]]@second$score) 
    return(df)
  })

# E11_5 methylation overlapped mESC shifting prom (mESC vs E16.5F) - 2442,  2513 out of 825 ranges
  overlap_E11_5F_methyl_signal.dfl <- lapply(1:length(overlap_E11_5F_WCG.l), function(x) {
    df <- data.frame(sample = names(overlap_E11_5F_WCG.l)[x], score = overlap_E11_5F_WCG.l[[x]]@second$score) 
    return(df)
  })

# E16_5F soma - methylation overlapped mESC shifting prom (mESC vs E16.5F) 1459, 1251, 1624  out of 825
  overlap_somaE16_5F_methyl_signal.dfl <- lapply(1:length(overlap_soma_E16_5F_WCG.l), function(x) {
    df <- data.frame(sample = names(overlap_soma_E16_5F_WCG.l)[x], score = overlap_soma_E16_5F_WCG.l[[x]]@second$score) 
    return(df)
  })

# E13_5F soma - methylation overlapped mESC shifting prom (mESC vs E16.5F) 1796, 2435, 2330  out of 825
  overlap_E13_5F_methyl_signal.dfl <- lapply(1:length(overlap_E13_5F_WCG.l), function(x) {
    df <- data.frame(sample = names(overlap_E13_5F_WCG.l)[x], score = overlap_E13_5F_WCG.l[[x]]@second$score) 
    return(df)
  })

# E16_5M  - methylation overlapped mESC shifting prom (mESC vs E16.5F) 1796, 2435, 2330  out of 825
  overlap_E16_5M_methyl_signal.dfl <- lapply(1:length(overlap_E16_5M_WCG.l), function(x) {
    df <- data.frame(sample = names(overlap_E16_5M_WCG.l)[x], score = overlap_E16_5M_WCG.l[[x]]@second$score) 
    return(df)
  })

# E13_5M - methylation overlapped mESC shifting prom (mESC vs E16.5F) 1796, 2435, 2330  out of 825
  overlap_E13_5M_methyl_signal.dfl <- lapply(1:length(overlap_E13_5M_WCG.l), function(x) {
    df <- data.frame(sample = names(overlap_E13_5M_WCG.l)[x], score = overlap_E13_5M_WCG.l[[x]]@second$score) 
    return(df)
  })

# bind into one dataframe
  overlap_E16_5F_methyl_signal.df <- do.call(rbind, overlap_E16_5F_methyl_signal.dfl)
  overlap_E16_5F_methyl_signal.df$sample <- factor(overlap_E16_5F_methyl_signal.df$sample, levels = names(overlap_E16_5F_WCG.l))
  
  overlap_E11_5F_methyl_signal.df <- do.call(rbind, overlap_E11_5F_methyl_signal.dfl)
  overlap_E11_5F_methyl_signal.df$sample <- factor(overlap_E11_5F_methyl_signal.df$sample, levels = names(overlap_E11_5F_WCG.l))
  
  overlap_somaE16_5F_methyl_signal.df <- do.call(rbind, overlap_somaE16_5F_methyl_signal.dfl)
  overlap_somaE16_5F_methyl_signal.df$sample <- factor(overlap_somaE16_5F_methyl_signal.df$sample, levels = names(overlap_soma_E16_5F_WCG.l))
  
  overlap_E13_5F_methyl_signal.df <- do.call(rbind, overlap_E13_5F_methyl_signal.dfl)
  overlap_E13_5F_methyl_signal.df$sample <- factor(overlap_E13_5F_methyl_signal.df$sample, levels = names(overlap_E13_5F_WCG.l))
  
  overlap_E16_5M_methyl_signal.df <- do.call(rbind, overlap_E16_5M_methyl_signal.dfl)
  overlap_E16_5M_methyl_signal.df$sample <- factor(overlap_E16_5M_methyl_signal.df$sample, levels = names(overlap_E16_5M_WCG.l))
  
  overlap_E13_5M_methyl_signal.df <- do.call(rbind, overlap_E13_5M_methyl_signal.dfl)
  overlap_E13_5M_methyl_signal.df$sample <- factor(overlap_E13_5M_methyl_signal.df$sample, levels = names(overlap_E13_5M_WCG.l))

# bind into one dataframe
  overlap_methyl_signal_narrower.df <- rbind(overlap_E11_5F_methyl_signal.df, overlap_E13_5F_methyl_signal.df, overlap_E13_5M_methyl_signal.df, overlap_E16_5F_methyl_signal.df, overlap_E16_5M_methyl_signal.df, overlap_somaE16_5F_methyl_signal.df)

# change levels
  overlap_methyl_signal_narrower.df$sample <- factor(overlap_methyl_signal_narrower.df$sample, levels = c("PGCE11_5F_WCG_rep1", "PGCE11_5F_WCG_rep2", 
                                                                                                          "PGCE13_5F_WCG_rep1", "PGCE13_5F_WCG_rep2", "PGCE13_5F_WCG_rep3",
                                                                                                          "PGCE13_5M_WCG_rep1", "PGCE13_5M_WCG_rep2", "PGCE13_5M_WCG_rep3",
                                                                                                          "PGCE16_5F_WCG_rep1", "PGCE16_5F_WCG_rep2", "PGCE16_5F_WCG_rep3",
                                                                                                          "PGCE16_5M_WCG_rep1", "PGCE16_5M_WCG_rep2",
                                                                                                          "soma_E16_5F_WCG_rep1", "soma_E16_5F_WCG_rep2", "soma_E16_5F_WCG_rep3"))
## -- plotting methylation signal -- ##
# violin plot
  samples <- c("PGC E11.5 r1", "PGC E11.5 r2", 
               "PGC E13.5F r1", "PGC E13.5F r2", "PGC E13.5F r2", 
               "PGC E13.5M r1", "PGC E13.5M r2", "PGC E13.5M r2", 
               "PGC E16.5F r1", "PGC E16.5F r2", "PGC E16.5F r3", 
               "PGC E16.5M r1", "PGC E16.5M r2",
               "soma E16.5F r1", "soma E16.5F r2", "soma E16.5F r3")

  p <- ggplot(overlap_methyl_signal_narrower.df, aes(y = score, x = sample, fill = sample), alpha = 0.6) +
    geom_violin(trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75), size = 0.5) +
    scale_fill_manual(values = c("#D1E11CFF", "#D1E11CFF",
                                 "#8BD646FF", "#8BD646FF", "#8BD646FF",
                                 "#75D054FF", "#75D054FF", "#75D054FF",
                                 "#3FBC73FF", "#3FBC73FF", "#3FBC73FF", 
                                 "#31B57BFF", "#31B57BFF",
                                 "#31B57BFF" , "#31B57BFF" , "#31B57BFF")) +
    scale_x_discrete(labels = samples) +
    theme_light() +
    theme(text = element_text(size = 18),
          axis.text.x = element_text(size = 18, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18)) +
    labs( x = "Samples") +
    guides(fill = F) +
    coord_cartesian(ylim = c(-0.2, 0.2))

  
    pdf("Figures/CpG_methylation_E16_5F_spread_diff_narrower.pdf", width = 8, height = 6)
      print(p)
    dev.off() 