#/*==========================================================================#*/
#' ## Extended Data Figure 16D
#/*==========================================================================#*/
# alternative GTF expression levels - Taf4b, Taf7l2, Taf7l, Taf9b

# load tag cluster data
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

# select for the transcript - Taf4b, Taf7l, Taf7l2, Taf9b
  Taf4b.grl <- lapply(tc_PGCs_oocyte_embryo_anno_symbol.grl, function(x) {
    return(x[grep("Taf4b", x$gene_symbol)])
  })

  Taf7l.grl <- lapply(tc_PGCs_oocyte_embryo_anno_symbol.grl, function(x) {
    return(x[grep("Taf7l", x$gene_symbol)])
  })

  Taf7l2.grl <- lapply(tc_PGCs_oocyte_embryo_anno_symbol.grl, function(x) {
    return(x[grep("Taf7l2", x$gene_symbol)])
  })

  Taf9b.grl <- lapply(tc_PGCs_oocyte_embryo_anno_symbol.grl, function(x) {
    return(x[grep("Taf9b", x$gene_symbol)])
  })

# exclude Toocyte, 2-cell and 4-cell
  Taf4b.grl <- Taf4b.grl[-c(13:20)]
  Taf7l.grl <- Taf7l.grl[-c(13:20)]
  Taf7l2.grl <- Taf7l2.grl[-c(13:20)]
  Taf9b.grl <- Taf9b.grl[-c(13:20)]


# filter to include only promoter regions and 5'UTR
  Taf4b_prom.grl <- lapply(Taf4b.grl, function(x) {
    x <- x[x$annotation == "Promoter" | x$annotation == "5' UTR"]
    return(x)
  })

  Taf7l_prom.grl <- lapply(Taf7l.grl, function(x) {
    x <- x[x$annotation == "Promoter" | x$annotation == "5' UTR"]
    return(x)
  })

  Taf7l2_prom.grl <- lapply(Taf7l2.grl, function(x) {
    x <- x[x$annotation == "Promoter" | x$annotation == "5' UTR"]
    return(x)
  })

  Taf9b_prom.grl <- lapply(Taf9b.grl, function(x) {
    x <- x[x$annotation == "Promoter" | x$annotation == "5' UTR"]
    return(x)
  })

# remove empty elemnts from taf7l
  Taf7l_prom.grl <- Taf7l_prom.grl[-c(1:3)]
  Taf7l2_prom.grl <- Taf7l2_prom.grl[-c(1:3)]

# convert to dataframes - GTFs set
  Taf4b_prom.dfl <- lapply(1:length(Taf4b_prom.grl), function(x) {
    df <- data.frame(sample = names(Taf4b_prom.grl)[[x]] , 
                     tpm =  Taf4b_prom.grl[[x]]$tpm, 
                     gene_symbol = Taf4b_prom.grl[[x]]$gene_symbol,
                     annotation = Taf4b_prom.grl[[x]]$annotation)
  })

  Taf7l_prom.dfl <- lapply(1:length(Taf7l_prom.grl), function(x) {
    df <- data.frame(sample = names(Taf7l_prom.grl)[[x]] , 
                     tpm =  Taf7l_prom.grl[[x]]$tpm, 
                     gene_symbol = Taf7l_prom.grl[[x]]$gene_symbol,
                     annotation = Taf7l_prom.grl[[x]]$annotation)
  })

  Taf7l2_prom.dfl <- lapply(1:length(Taf7l2_prom.grl), function(x) {
    df <- data.frame(sample = names(Taf7l2_prom.grl)[[x]] , 
                     tpm =  Taf7l2_prom.grl[[x]]$tpm, 
                     gene_symbol = Taf7l2_prom.grl[[x]]$gene_symbol,
                     annotation = Taf7l2_prom.grl[[x]]$annotation)
  }) 

  Taf9b_prom.dfl <- lapply(1:length(Taf9b_prom.grl), function(x) {
    df <- data.frame(sample = names(Taf9b_prom.grl)[[x]] , 
                     tpm =  Taf9b_prom.grl[[x]]$tpm, 
                     gene_symbol = Taf9b_prom.grl[[x]]$gene_symbol,
                     annotation = Taf9b_prom.grl[[x]]$annotation)
  })     

  names(Taf4b_prom.dfl) <- names(Taf4b_prom.grl)
  names(Taf7l_prom.dfl) <- names(Taf7l_prom.grl)
  names(Taf7l2_prom.dfl) <- names(Taf7l2_prom.grl)
  names(Taf9b_prom.dfl) <- names(Taf9b_prom.grl)


# sum tpms of same genes
  library(magrittr)
  library(dplyr)
  
  Taf4b_prom_summed.dfl <- lapply(1:length(Taf4b_prom.dfl), function(x) {
    df <- Taf4b_prom.dfl[[x]] %>% dplyr::group_by(gene_symbol) %>% dplyr::select(tpm) %>% dplyr::summarise(sum_tpm = sum(tpm))
    df <- as.data.frame(df)
    df$sample <- names(Taf4b_prom.dfl)[[x]]
    return(df)
  })
  names(Taf4b_prom_summed.dfl) <- names(Taf4b_prom.dfl)

  Taf7l_prom_summed.dfl <- lapply(1:length(Taf7l_prom.dfl), function(x) {
    df <- Taf7l_prom.dfl[[x]] %>% dplyr::group_by(gene_symbol) %>% dplyr::select(tpm) %>% dplyr::summarise(sum_tpm = sum(tpm))
    df <- as.data.frame(df)
    df$sample <- names(Taf7l_prom.dfl)[[x]]
    return(df)
  })  
  names(Taf7l_prom_summed.dfl) <- names(Taf7l_prom.dfl)


  Taf7l2_prom_summed.dfl <- lapply(1:length(Taf7l2_prom.dfl), function(x) {
    df <- Taf7l2_prom.dfl[[x]] %>% dplyr::group_by(gene_symbol) %>% dplyr::select(tpm) %>% dplyr::summarise(sum_tpm = sum(tpm))
    df <- as.data.frame(df)
    df$sample <- names(Taf7l2_prom.dfl)[[x]]
    return(df)
  })  
  names(Taf7l2_prom_summed.dfl) <- names(Taf7l2_prom.dfl)

  Taf9b_prom_summed.dfl <- lapply(1:length(Taf9b_prom.dfl), function(x) {
    df <- Taf9b_prom.dfl[[x]] %>% dplyr::group_by(gene_symbol) %>% dplyr::select(tpm) %>% dplyr::summarise(sum_tpm = sum(tpm))
    df <- as.data.frame(df)
    df$sample <- names(Taf9b_prom.dfl)[[x]]
    return(df)
  })  
  names(Taf9b_prom_summed.dfl) <- names(Taf9b_prom.dfl)


# convert all lists to dataframes
  Taf4b_prom_summed.df <- do.call(rbind, Taf4b_prom_summed.dfl)
  Taf7l_prom_summed.df <- do.call(rbind, Taf7l_prom_summed.dfl)
  Taf7l2_prom_summed.df <- do.call(rbind, Taf7l2_prom_summed.dfl)
  Taf9b_prom_summed.df <- do.call(rbind, Taf9b_prom_summed.dfl)


# combine all genes
  GTF_all_summed.df <- rbind(Taf4b_prom_summed.df, Taf7l_prom_summed.df, Taf7l2_prom_summed.df, Taf9b_prom_summed.df)

# check levels
  GTF_all_summed.df$gene_symbol <- droplevels(GTF_all_summed.df$gene_symbol)
  GTF_all_summed.df$sample <- factor(GTF_all_summed.df$sample, levels = names(tc_PGCs_oocyte_embryo_anno.grl)[1:12])

# subselect only all and female
  library(dplyr)
  GTF_all_summed_fem.df <- GTF_all_summed.df %>% filter(sample %in% c("E14_mESC", "E9_5", "E10_5", "E11_5", "E12_5F", "E13_5F", "E14_5F", "E16_5F"))
  GTF_all_summed_fem.df$sample <- droplevels.factor(GTF_all_summed_fem.df$sample)

# subselect only all and male
  library(dplyr)
  GTF_all_summed_m.df <- GTF_all_summed.df %>% filter(sample %in% c("E14_mESC", "E9_5", "E10_5", "E11_5", "E12_5M", "E13_5M", "E14_5M", "E16_5M"))
  GTF_all_summed_m.df$sample <- droplevels.factor(GTF_all_summed_m.df$sample)    


library(ggplot2)
  pF <- ggplot(data = GTF_all_summed_fem.df, aes(x = sample,
                                                 y = sum_tpm,
                                                 fill = gene_symbol,
                                                 color = gene_symbol,
                                                 group = gene_symbol)) + 
    geom_line() +
    scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")) +
    geom_point(size = 4, shape = 21, colour = "black") + 
    scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")) +
    theme_light() +
    theme(text = element_text(size = 12, colour = "black"), 
          axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 12, colour = "black"),
          legend.position = "right") +
    scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, by = 50))
  
  pM <- ggplot(data = GTF_all_summed_m.df, aes(x = sample,
                                               y = sum_tpm,
                                               fill = gene_symbol,
                                               color = gene_symbol,
                                               group = gene_symbol)) + 
    geom_line() +
    scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")) +
    geom_point(size = 4, shape = 21, colour = "black") + 
    scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")) +
    theme_light() +
    theme(text = element_text(size = 12, colour = "black"), 
          axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 12, colour = "black"),
          legend.position = "right") +
    scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, by = 50))
  
  # combine figures
    library(cowplot)
    plot_grid(pF, pM, align = "h", ncol = 2)
    
  
    pdf("Figures/GTF_expression.pdf", width = 12, height = 4)
      plot_grid(pF, pM, align = "h", ncol = 2)
    dev.off() 