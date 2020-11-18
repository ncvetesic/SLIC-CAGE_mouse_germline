#/*==========================================================================#*/
#' ## Extended Data Figure 16C
#/*==========================================================================#*/
# H3.3 KO analysis - E16.5 shifts in the context of differentially expressed genes upon H3.3 KO

# read in table of diff expressed genes (H3.3 KO vs WT)
  library(readxl)
  diff_exp_H3_3.df <- read_xls("/Users/ncvetesic/Documents/projects/PGCs_mouse/extra_data/H3_3/supp_29.13.1377_TableS2.xls", skip = 1, sheet = 1, col_types = c("numeric", "numeric", "numeric", "text", "text", "text", "text", "text", "text", "text"))

# separate upreg and downreg
  diff_exp_H3_3_upreg.df <- diff_exp_H3_3.df[diff_exp_H3_3.df$`Log Ratio` >= 0.5, ]
  diff_exp_H3_3_downreg.df <- diff_exp_H3_3.df[diff_exp_H3_3.df$`Log Ratio` <= -0.5, ]

# filter for p.value
  diff_exp_H3_3_upreg_filtr.df <- diff_exp_H3_3_upreg.df[diff_exp_H3_3_upreg.df$`False Discovery Rate (q-value)` <= 0.05, ]
  diff_exp_H3_3_downreg_filtr.df <- diff_exp_H3_3_downreg.df[diff_exp_H3_3_downreg.df$`False Discovery Rate (q-value)` <= 0.05, ]

# import shitfing promoters annotated
  PGC_oocyte_embryo_vs_mESC_shift_anno.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/PGC_oocyte_embryo_vs_mESC_shift_anno_grl.RDS")

# convert entrez id's to gene symbols
  library(org.Mm.eg.db)
  library(AnnotationDbi)

# extract ENTREZ ids from annotated tag clusters 
  PGC_oocyte_embryo_vs_mESC_shift_anno.grl <- lapply(PGC_oocyte_embryo_vs_mESC_shift_anno.grl, function(x) {
    gene_ids <- mcols(x)$geneId
  
  # get all gene symbols from org.Mm.eg.db
    geneSymbols <- mapIds(org.Mm.eg.db, keys = gene_ids, column="SYMBOL", keytype="ENTREZID", multiVals= "first")
    geneEnsembl <- mapIds(org.Mm.eg.db, keys = gene_ids, column="ENSEMBL", keytype="ENTREZID", multiVals= "first")
  
  # attach gene symbols to annotated tag clusters
    x$ENSEMBL <- as.character(geneEnsembl)
    x$gene_symbol <- as.character(geneSymbols)  
    return(x)
  })
  sum(PGC_oocyte_embryo_vs_mESC_shift_anno.grl[["E16_5F"]]$gene_symbol %in% diff_exp_H3_3_upreg_filtr.df$ID) # 7 overlapping

  sum(PGC_oocyte_embryo_vs_mESC_shift_anno.grl[["E16_5F"]]$gene_symbol %in% diff_exp_H3_3_downreg_filtr.df$ID) # 5 overlaping

# check with unfiltered set not many may be expressed...
  sum(PGC_oocyte_embryo_vs_mESC_shift_anno.grl[["E16_5F"]]$gene_symbol %in% diff_exp_H3_3.df$ID)

# ---- select all overlaping genes -  747 is found - check log2 distribution of overlapping
  diff_exp_H3_3.df$`Log Ratio`[diff_exp_H3_3.df$ID %in% PGC_oocyte_embryo_vs_mESC_shift_anno.grl[["E16_5F"]]$gene_symbol]
  hist(diff_exp_H3_3.df$`Log Ratio`[diff_exp_H3_3.df$ID %in% PGC_oocyte_embryo_vs_mESC_shift_anno.grl[["E16_5F"]]$gene_symbol], 200)

# overlap distributions of overlap and all diff exp genes
  all.df <- data.frame("sample" = "all", "log2_ratio" = diff_exp_H3_3.df$`Log Ratio`, "q_value" = diff_exp_H3_3.df$`False Discovery Rate (q-value)`)
  table_all.df <- data.frame("percentage" = table(all.df$log2_ratio)/sum(table(all.df$log2_ratio))*100)
  table_all.df$sample <- "all"
  colnames(table_all.df) <- c("log2_ratio", "percentage", "sample")
  table_all.df$log2_ratio <- as.numeric(names(table(all.df$log2_ratio)))

  overlap.df <- data.frame("sample" = "overlap", "log2_ratio" = diff_exp_H3_3.df$`Log Ratio`[diff_exp_H3_3.df$ID %in% PGC_oocyte_embryo_vs_mESC_shift_anno.grl[["E16_5F"]]$gene_symbol], "q_value" = diff_exp_H3_3.df$`False Discovery Rate (q-value)`[diff_exp_H3_3.df$ID %in% PGC_oocyte_embryo_vs_mESC_shift_anno.grl[["E16_5F"]]$gene_symbol])
  table_overlap.df <- data.frame("percentage" = table(overlap.df$log2_ratio)/sum(table(overlap.df$log2_ratio))*100)
  table_overlap.df$sample <- "overlap"
  colnames(table_overlap.df) <- c("log2_ratio", "percentage", "sample")
  table_overlap.df$log2_ratio <- as.numeric(names(table(overlap.df$log2_ratio)))

# combine both dataframes
  table_both.df <- rbind(table_all.df, table_overlap.df)

# change levels
  table_both.df$sample <- factor(table_both.df$sample, levels = c("all", "overlap"))
  table_both.df$log2_ratio <- as.numeric(table_both.df$log2_ratio)

# plot distribution 
  library(ggplot2)
  p <- ggplot(table_both.df, aes(x = log2_ratio, fill = sample, alpha = 0.85)) +
    geom_bar(col = "black", stat = "bin", lwd = 0.5, alpha = 0.4, binwidth = 0.1) +
    scale_fill_manual("sample", values = c("firebrick1", "#3FBC73FF")) +
    theme_light()  +
    theme(text = element_text(size = 16, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          legend.position = "top") +
    labs(y = "%", x = "log2 ratio") +
    scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 0.5)) 


# calculate median of overlaping genes and randomly selected genes:
  med_overlap <- median(overlap.df$log2_ratio)

# select - 741 random diff exp genes
  med_random <- vector()
  set.seed(7)
  
  for(i in 1:10000) {
    log2_random <- all.df$log2_ratio[sample(1:nrow(all.df), size = 741)]
    med_random[i] <- median(log2_random)
  }

# create dataframe
  medians.df <- data.frame("medians" = med_random)  
  
  p <- ggplot(medians.df, aes(x = medians, alpha = 0.85)) +
    geom_bar(col = "black", stat = "bin", lwd = 0.5, alpha = 0.4, binwidth = 0.001) +
    scale_fill_manual(values = c("#3FBC73FF")) +
    theme_light()  +
    theme(text = element_text(size = 16, colour = "black"),
          legend.title = element_blank(),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          legend.position = "top") +
    labs(y = "Counts", x = "log2 ratio") +
    geom_vline(xintercept = median(medians.df$medians), col = "black", linetype = "dashed", size = 1) +
    geom_vline(xintercept = med_overlap, col = "firebrick", linetype = "dashed", size = 1) +
    scale_x_continuous(limits = c(-0.075, 0.075), breaks = seq(-0.075, 0.075, by = 0.05)) +
    scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, by = 200)) 


  pdf("Figures/H3_3logratio_diffexprs_significance.pdf", height = 4, width = 8)
    print(p)
  dev.off()   