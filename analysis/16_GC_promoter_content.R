#/*==========================================================================#*/
#' ## Extended Data Figure 4C, 12C
#/*==========================================================================#*/
# GC content of shifting promoters 

# import data
  library(BSgenome.Mmusculus.UCSC.mm10)
  domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/domTSS_PGC_oocyte_embryo_vs_mESC_shift_grl.RDS")
  domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/domTSS_PGC_oocyte_embryo_vs_mESC_shift_shift_X_grl.RDS")
  PGC_oocyte_embryo_vs_mESC_shift_anno.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/PGC_oocyte_embryo_vs_mESC_shift_anno_grl.RDS")

  domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl <- lapply(domTSS_PGC_oocyte_embryo_vs_mESC_shift.grl, function(x){
    gr <- dropSeqlevels(x, "chrM", pruning.mode = "coarse")
  })

  domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl <- lapply(domTSS_PGC_oocyte_embryo_vs_mESC_shift_X.grl, function(x){
    gr <- dropSeqlevels(x, "chrM", pruning.mode = "coarse")
  })

  PGC_oocyte_embryo_vs_mESC_shift_anno.grl <- lapply(PGC_oocyte_embryo_vs_mESC_shift_anno.grl, function(x){
    gr <- dropSeqlevels(x, "chrM", pruning.mode = "coarse")
  })

# load annotated all TCs - I need cluster tpm info 
  tc_PGCs_oocyte_embryo_anno.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/tc_PGCs_oocyte_embryo_anno_grl.RDS")

  tc_PGCs_oocyte_embryo_anno.grl <- lapply(tc_PGCs_oocyte_embryo_anno.grl, function(x){
    gr <- dropSeqlevels(x, "chrM", pruning.mode = "coarse")
  })

# use just the shifting consensus cluster
# extract sequence
  shifting_proms_seq.l <- lapply(PGC_oocyte_embryo_vs_mESC_shift_anno.grl, function(x) {
    return(getSeq(BSgenome.Mmusculus.UCSC.mm10, x))
  })

# count CpG in each and divide by width
  shifting_CpG_count.l <- lapply(1:length(shifting_proms_seq.l), function(x) {
    counts <- dinucleotideFrequency(shifting_proms_seq.l[[x]], step = 1, as.prob = F)
    CpG_count.df <- data.frame(CG_GC_counts = counts[, "CG"] + counts[, "GC"], 
                               width = width(shifting_proms_seq.l[[x]]), 
                               norm_CG_GC = (counts[, "CG"] + counts[, "GC"])/(width(shifting_proms_seq.l[[x]])-1),
                               sample = names(PGC_oocyte_embryo_vs_mESC_shift_anno.grl)[x])
    CpG_count.df <- CpG_count.df[!is.na(CpG_count.df$norm_CG_GC), ]
    
    return(CpG_count.df)
  })

# calculate CG_GC content in consensus promoters in general
  # extract consensus promoters
    mouse_all_consensus.df <- consensusClusters(CAGEset_PGC_embryo_merged)

# turn consensus clusters into granges
  mouse_all_consensus.gr <- GRanges(seqnames = mouse_all_consensus.df$chr,
                                    ranges = IRanges(start = mouse_all_consensus.df$start,
                                                     end = mouse_all_consensus.df$end),
                                    strand = mouse_all_consensus.df$strand,
                                    tpm = mouse_all_consensus.df$tpm)

seqinfo(mouse_all_consensus.gr) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(mouse_all_consensus.gr)]

# save consensus cluster GRanges
  saveRDS(mouse_all_consensus.gr, "merged_replicates/intermediate_data/mouse_all_consensus_gr.RDS")

# consensus cluster cpg count
  # extract sequence
    consensus_counts_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, mouse_all_consensus.gr)
    counts <- dinucleotideFrequency(consensus_counts_seq, step = 1, as.prob = F)
    CpG_consensus_CG_GC_counts.df <- data.frame(CG_GC_counts = counts[, "CG"] + counts[, "GC"],
                                                width = width(consensus_counts_seq),
                                                norm_CG_GC = (counts[, "CG"] + counts[, "GC"])/(width(consensus_counts_seq)-1))
    CpG_consensus_CG_GC_counts.df <- CpG_consensus_CG_GC_counts.df[!is.na(CpG_consensus_CG_GC_counts.df$norm_CG_GC), ]
    CpG_consensus_CG_GC_counts.df$sample <- "consensus"

# combine both dataframes - consensus and sample
  all_CpG_count.l <- c(shifting_CpG_count.l, list(CpG_consensus_CG_GC_counts.df))  

# flatten to a dataframe
  all_CpG_count.df <- do.call(rbind, all_CpG_count.l)

# ------ select 16.5F/M, PN6, oocyte, 2-cell, 4-cell and consensus
  all_CpG_count_sel.df <- all_CpG_count.df[all_CpG_count.df$sample == "E16_5F" | all_CpG_count.df$sample == "E16_5M" | all_CpG_count.df$sample == "PN6" | all_CpG_count.df$sample == "PN14" | all_CpG_count.df$sample == "oocyte" | all_CpG_count.df$sample == "S2_cell" | all_CpG_count.df$sample == "S4_cell" | all_CpG_count.df$sample == "consensus", ]

# drop unused sequence levels
  all_CpG_count_sel.df$sample <- droplevels(all_CpG_count_sel.df$sample)

col <- col_scheme[c("E16.5F", "E16.5M", "P6", "P14", "MII", "2-cell", "4-cell")]
samples <- levels(all_CpG_count_sel.df$sample)

# ggplot 
  p <- ggplot(all_CpG_count_sel.df, aes(y = norm_CG_GC, x = sample , fill = sample, alpha = 0.6)) +
    geom_boxplot(outlier.shape = NA, varwidth = T, notch = F)  +
    scale_fill_manual(values = c("#3FBC73FF", "#31B57BFF", "#440154FF", "#46337EFF", "#3F4889FF", "#365C8DFF", "#2E6E8EFF", "gray66")) +
    scale_x_discrete(labels = samples) +
    theme_bw() +
    theme(text = element_text(size = 16, color = "black"),
          axis.text.x = element_text(size = 16, angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(size = 16, color = "black"),
          axis.title.x = element_text(size = 16, color = "black"),
          axis.title.y = element_text(size = 16,color = "black"),
          legend.position = "none") +
    labs(title = "") +
    labs( x = "", y = "CG norm") +
    guides(fill = FALSE) +
    scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2)) +
    geom_hline(yintercept = median(all_CpG_count_sel.df$norm_CG_GC), linetype="dashed", color = "indianred1", size=1, alpha = 0.8)

  pdf("Figures/shifting_cpg_boxplots.pdf", height = 4, width = 4)
    print(p)
  dev.off()
  
  # ----- GC content in expanded shifting promoters ----- # 
  # expand consensus clusters to include 1000 bp upstream and downstream
    PGC_oocyte_embryo_vs_mESC_shift_anno_expanded.grl <- lapply(PGC_oocyte_embryo_vs_mESC_shift_anno.grl, function(x) {
      return(resize(x, 2000, fix = "center", use.names = T, ignore.strand = FALSE))
    })
  
  # extract sequence
    shifting_proms_seq.l <- lapply(PGC_oocyte_embryo_vs_mESC_shift_anno_expanded.grl, function(x) {
      return(getSeq(BSgenome.Mmusculus.UCSC.mm10, x))
    })
  
  # count CpG in each and divide by width
    shifting_CpG_count.l <- lapply(1:length(shifting_proms_seq.l), function(x) {
      counts <- dinucleotideFrequency(shifting_proms_seq.l[[x]], step = 1, as.prob = F)
      CpG_count.df <- data.frame(CG_GC_counts = counts[, "CG"] + counts[, "GC"], 
                                 width = width(shifting_proms_seq.l[[x]]), 
                                 norm_CG_GC = (counts[, "CG"] + counts[, "GC"])/(width(shifting_proms_seq.l[[x]])-1),
                                 sample = names(PGC_oocyte_embryo_vs_mESC_shift_anno.grl)[x])
      CpG_count.df <- CpG_count.df[!is.na(CpG_count.df$norm_CG_GC), ]
      
      return(CpG_count.df)
    })
  
  
  # expand consensus promoters
    mouse_all_consensus_expanded.gr <- resize(mouse_all_consensus.gr, 2000, fix = "center", use.names = T, ignore.strand = FALSE)
    mouse_all_consensus_expanded.gr <- mouse_all_consensus_expanded.gr[width(trim(mouse_all_consensus_expanded.gr)) == width(mouse_all_consensus_expanded.gr)]
  
  # consensus cluster cpg count
    # extract sequence
      consensus_counts_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, mouse_all_consensus_expanded.gr)
      counts <- dinucleotideFrequency(consensus_counts_seq, step = 1, as.prob = F)
      CpG_consensus_CG_GC_counts.df <- data.frame(CG_GC_counts = counts[, "CG"] + counts[, "GC"],
                                                  width = width(consensus_counts_seq),
                                                  norm_CG_GC = (counts[, "CG"] + counts[, "GC"])/(width(consensus_counts_seq)-1))
      CpG_consensus_CG_GC_counts.df <- CpG_consensus_CG_GC_counts.df[!is.na(CpG_consensus_CG_GC_counts.df$norm_CG_GC), ]
      CpG_consensus_CG_GC_counts.df$sample <- "consensus"
  
  # combine both dataframes - consensus and sample
    all_CpG_count.l <- c(shifting_CpG_count.l, list(CpG_consensus_CG_GC_counts.df))  
  
  # flatten to a dataframe
    all_CpG_count.df <- do.call(rbind, all_CpG_count.l)
  
  # ---- select 16.5F/M, PN6, PN7, oocyte, 2-cell, 4-cell and consensus
    all_CpG_count_sel.df <- all_CpG_count.df[all_CpG_count.df$sample == "E16_5F" | all_CpG_count.df$sample == "E16_5M" | all_CpG_count.df$sample == "PN6" | all_CpG_count.df$sample == "PN14" | all_CpG_count.df$sample == "oocyte" | all_CpG_count.df$sample == "S2_cell" | all_CpG_count.df$sample == "S4_cell" | all_CpG_count.df$sample == "consensus", ]
  
  # drop unused sequence levels
    all_CpG_count_sel.df$sample <- droplevels(all_CpG_count_sel.df$sample)
  
  col <- col_scheme[c("E16.5F", "E16.5M", "P6", "P14", "MII", "2-cell", "4-cell")]
  samples <- levels(all_CpG_count_sel.df$sample)
  
  # ggplot
    p <- ggplot(all_CpG_count_sel.df, aes(y = norm_CG_GC, x = sample , fill = sample, alpha = 0.6)) +
      geom_boxplot(outlier.shape = NA, varwidth = T, notch = F)  +
      scale_fill_manual(values = c("#3FBC73FF", "#31B57BFF", "#440154FF", "#46337EFF", "#3F4889FF", "#365C8DFF", "#2E6E8EFF", "gray66")) +
      scale_x_discrete(labels = samples) +
      theme_bw() +
      theme(text = element_text(size = 16, color = "black"),
            axis.text.x = element_text(size = 16, angle = 45, hjust = 1, color = "black"),
            axis.text.y = element_text(size = 16, color = "black"),
            axis.title.x = element_text(size = 16, color = "black"),
            axis.title.y = element_text(size = 16,color = "black"),
            legend.position = "none") +
      labs(title = "") +
      labs( x = "", y = "CG norm") +
      guides(fill = FALSE) +
      scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
      geom_hline(yintercept = median(all_CpG_count_sel.df$norm_CG_GC), linetype="dashed", color = "indianred1", size=1, alpha = 0.8)
  
    pdf("Figures/shifting_expanded_cpg_boxplots.pdf", height = 4, width = 4)
      print(p)
    dev.off()
  
