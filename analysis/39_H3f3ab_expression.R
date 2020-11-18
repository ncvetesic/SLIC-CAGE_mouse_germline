#/*==========================================================================#*/
#' ## Extended Data Figure 16B
#/*==========================================================================#*/
# H3f3a and H3f3b expression

# import tag cluster data
  tc_PGCs_oocyte_embryo_anno.grl <- readRDS( "../all_reps_PN6_2cell/intermediate_data/tc_PGCs_oocyte_embryo_anno_grl.RDS")

# add gene symbol annotations - so I can grep for hist
# convert entrez id's to gene symbols
  library(org.Mm.eg.db)
  library(AnnotationDbi)

# extract ENTREZ ids from annotated tag clusters 
  tc_PGCs_oocyte_embryo_anno_symbol.grl <- lapply(tc_PGCs_oocyte_embryo_anno.grl, function(x) {
    gene_ids <- mcols(x)$geneId
    
    # get all gene symbols from org.Mm.eg.db
    geneSymbols <- mapIds(org.Mm.eg.db, keys = gene_ids, column="SYMBOL", keytype="ENTREZID", multiVals= "first")
    geneEnsembl <- mapIds(org.Mm.eg.db, keys = gene_ids, column="ENSEMBL", keytype="ENTREZID", multiVals= "first")
    
    # attach gene symbols to annotated tag clusters
    x$ENSEMBL <- as.character(geneEnsembl)
    x$gene_symbol <- as.character(geneSymbols)  
    return(x)
  })

# select for the transcript
  histone.grl <- lapply(tc_PGCs_oocyte_embryo_anno_symbol.grl, function(x) {
    return(x[grep("Hist", x$gene_symbol)])
  })

# exclude TBP2KO and PN7
  histone.grl <- histone.grl[-c(14, 15, 17)]

# filter to include only promoter regions and 5'UTR
  histone_prom.grl <- lapply(histone.grl, function(x) {
    x <- x[x$annotation == "Promoter" | x$annotation == "5' UTR"]
    return(x)
  })

# select H3f3a and H3f3b
  H3f3ab.grl <- lapply(tc_PGCs_oocyte_embryo_anno_symbol.grl, function(x) {
    x <- x[grep("H3f3", x$gene_symbol)]
    return(x)
  })

# filter to include only promoter regions and 5'UTR
  H3f3ab_prom.grl <- lapply(H3f3ab.grl, function(x) {
    x <- x[x$annotation == "Promoter" | x$annotation == "5' UTR"]
    return(x)
  })

# convert to dataframes - histone set
  histone_prom.dfl <- lapply(1:length(histone_prom.grl), function(x) {
    df <- data.frame(sample = names(histone_prom.grl)[[x]] , 
                     tpm =  histone_prom.grl[[x]]$tpm, 
                     gene_symbol = histone_prom.grl[[x]]$gene_symbol,
                     annotation = histone_prom.grl[[x]]$annotation)
  })

# convert to dataframes - h3f3ab
  H3f3ab_prom.dfl <- lapply(1:length(H3f3ab_prom.grl), function(x) {
    df <- data.frame(sample = names(H3f3ab_prom.grl)[[x]] , 
                     tpm =  H3f3ab_prom.grl[[x]]$tpm, 
                     gene_symbol = H3f3ab_prom.grl[[x]]$gene_symbol,
                     annotation = H3f3ab_prom.grl[[x]]$annotation)
  })
  
  names(H3f3ab_prom.dfl) <- names(H3f3ab_prom.grl)

# exclude P7, P7 TBP2KO and P14 TBP2KO
  H3f3ab_prom.dfl <- H3f3ab_prom.dfl[-c(14, 15, 17)]

# sum tpms of same genes
  library(magrittr)
  library(dplyr)
  names(histone_prom.dfl) <- names(histone.grl)

  histone_prom_summed.dfl <- lapply(1:length(histone_prom.dfl), function(x) {
    df <- histone_prom.dfl[[x]] %>% dplyr::group_by(gene_symbol) %>% dplyr::select(tpm) %>% dplyr::summarise(sum_tpm = sum(tpm))
    df <- as.data.frame(df)
    df$sample <- names(histone_prom.dfl)[[x]]
    return(df)
  })

  H3f3ab_prom_summed.dfl <- lapply(1:length(H3f3ab_prom.dfl), function(x) {
    df <- H3f3ab_prom.dfl[[x]] %>% dplyr::group_by(gene_symbol) %>% dplyr::select(tpm) %>% dplyr::summarise(sum_tpm = sum(tpm))
    df <- as.data.frame(df)
    df$sample <- names(H3f3ab_prom.dfl)[[x]]
    return(df)
  })

# convert all lists to dataframes
  histone_prom_summed.df <- do.call(rbind, histone_prom_summed.dfl)
  H3f3ab_prom_summed.df <- do.call(rbind, H3f3ab_prom_summed.dfl)
  
  all_histone_prom_summed.df <- rbind(histone_prom_summed.df, H3f3ab_prom_summed.df)

# select histone3 variants
  hist3.df <- rbind(all_histone_prom_summed.df[grep("Hist3", all_histone_prom_summed.df$gene_symbol), ], all_histone_prom_summed.df[grep("H3f3", all_histone_prom_summed.df$gene_symbol), ])

# check levels
  hist3.df$gene_symbol <- droplevels(hist3.df$gene_symbol)
  hist3.df$sample <- factor(hist3.df$sample, levels = names(H3f3ab_prom.dfl))

  library(ggplot2)
  p <- ggplot(data = hist3.df, aes(x = sample,
                                   y = log2(sum_tpm + 1),
                                   fill = gene_symbol,
                                   color = gene_symbol,
                                   group = gene_symbol)) + 
    geom_line() +
    scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb","#e78ac3","#a6d854")) +
    geom_point(size = 4, shape = 21, colour = "black") + 
    scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb","#e78ac3","#a6d854")) +
    theme_light() +
    theme(text = element_text(size = 12, colour = "black"), 
          axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 12, colour = "black"),
          legend.position = "right") 
  
  
  pdf("Figures/hist3_expression.pdf", width = 8, height = 6)
    print(p)
  dev.off() 

# ------- print only h3h3 variants for SUpplementary Figure14
# check levels
  H3f3ab_prom_summed.df$sample <- factor(H3f3ab_prom_summed.df$sample, levels = names(H3f3ab_prom.dfl))

# subselect only all and female
  library(dplyr)
  H3f3ab_prom_summed_fem.df <- H3f3ab_prom_summed.df %>% filter(sample %in% c("E14_mESC", "E9_5", "E10_5", "E11_5", "E12_5F", "E13_5F", "E14_5F", "E16_5F"))
  H3f3ab_prom_summed_fem.df$sample <- droplevels.factor(H3f3ab_prom_summed_fem.df$sample)

# subselect only all and male
  library(dplyr)
  H3f3ab_prom_summed_m.df <- H3f3ab_prom_summed.df %>% filter(sample %in% c("E14_mESC", "E9_5", "E10_5", "E11_5", "E12_5M", "E13_5M", "E14_5M", "E16_5M"))
  H3f3ab_prom_summed_m.df$sample <- droplevels.factor(H3f3ab_prom_summed_m.df$sample)    

  library(ggplot2)
  pF <- ggplot(data = H3f3ab_prom_summed_fem.df, aes(x = sample,
                                                     y = sum_tpm,
                                                     fill = gene_symbol,
                                                     color = gene_symbol,
                                                     group = gene_symbol)) + 
    geom_line() +
    scale_color_manual(values = c("#fc8d62", "#8da0cb")) +
    geom_point(size = 4, shape = 21, colour = "black") + 
    scale_fill_manual(values = c("#fc8d62", "#8da0cb")) +
    theme_light() +
    theme(text = element_text(size = 12, colour = "black"), 
          axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 12, colour = "black"),
          legend.position = "right") +
    scale_y_continuous(limits = c(0, 1500), breaks = seq(0, 1500, by =500))
  
  pM <- ggplot(data = H3f3ab_prom_summed_m.df, aes(x = sample,
                                                   y = sum_tpm,
                                                   fill = gene_symbol,
                                                   color = gene_symbol,
                                                   group = gene_symbol)) + 
    geom_line() +
    scale_color_manual(values = c("#fc8d62", "#8da0cb")) +
    geom_point(size = 4, shape = 21, colour = "black") + 
    scale_fill_manual(values = c("#fc8d62", "#8da0cb")) +
    theme_light() +
    theme(text = element_text(size = 12, colour = "black"), 
          axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 12, colour = "black"),
          legend.position = "right") +
    scale_y_continuous(limits = c(0, 1500), breaks = seq(0, 1500, by =500))

# combine figures
  library(cowplot)
  plot_grid(pF, pM, align = "h", ncol = 2)


  pdf("Figures/hist3_3_expression.pdf", width = 12, height = 4)
    plot_grid(pF, pM, align = "h", ncol = 2)
  dev.off() 