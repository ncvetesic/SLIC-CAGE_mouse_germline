#/*==========================================================================#*/
#' ## Extended Data Figure 2E,F
#/*==========================================================================#*/
# PCA plot at CTSS level using DeSEQ2

# need the raw TPM values, so I need to set the normalize method to "none", recreate consensus clusters and extract TPM - see 02_TC_tpm_PCA.R

# load raw CAGE set
  rawCAGEset_merged_reps <- readRDS(file = "../all_reps_PN6_2cell/intermediate_data/rawCAGEset_merged_reps.RDS")

## Prepare the raw count table
# extract the raw counts
  CTSS_raw.counts.df <- rawCAGEset_merged_reps@normalizedTpmMatrix

# prepare DeSeq2 object
library(DESeq2)

# create a count matrix for DESeq2
count.matrix <- as.matrix(CTSS_raw.counts.df)

# info table/ design table / colData
  samples <- colnames(count.matrix)
  info.df <- data.frame(condition = factor(x = c(rep("PGC", times = 11), rep("oocyte", times  = 4), rep("embryo", times  = 2)), 
                                           levels = c("PGC", "oocyte", "embryo")),
                        row.names = samples)

# create DESeqDataSet object - stores count matrix and design
  dds <- DESeqDataSetFromMatrix(countData = count.matrix,
                                colData = info.df,
                                design = ~ condition)

# PCA preparation
  dds <- DESeq(dds)
  vsd <- vst(dds, blind = TRUE) # blind means don't take design into consideration

# save dds object PGC embryo oocyte  
  saveRDS (dds, "../all_reps_PN6_2cell/intermediate_data/CTSS_dds_all_merg_no_mESC_no_P7_noTBP2.RDS")

# using top 8,000,000 most variable CTSSs (out of 10 mil)
  pcaData <- plotPCA(vsd, intgroup=c("condition"), ntop = 8000000, returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

# plotting of PCA
  library(ggrepel)

# change names of samples and levels
  samples <- factor(c("mESC-E14", "E9.5", "E10.5", "E11.5", "E12.5F", "E12.5M", "E13.5F", "E13.5M", "E14.5F", "E14.5M", "E16.5F", "E16.5M", "P6", "P14", "MII", "2-cell", "4-cell"), levels = c("mESC-E14", "E9.5", "E10.5", "E11.5", "E12.5F", "E12.5M", "E13.5F", "E13.5M", "E14.5F", "E14.5M", "E16.5F", "E16.5M", "P6", "P14", "MII", "2-cell",  "4-cell"))
  pcaData$name <- samples

# setting plotting colour
  col <- col_scheme

# plotting 
  p <-ggplot(pcaData, aes(PC1, PC2, fill = name)) +
    geom_point(size=3, pch = 21, colour = "black") +
    scale_fill_manual(values = col) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    theme_bw() +
    theme(text = element_text(size = 20, colour = "black"), 
          axis.text.x = element_text(size = 18, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black"),
          legend.title = element_blank(),
          legend.position = "none") + 
    geom_label_repel(aes(label = name),
                     fill = "white", colour = "grey50",
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50',
                     show.legend = F) +
    scale_x_continuous(limits = c(-1200, 2200), breaks = seq(-1200, 2200, by = 1000)) +
    scale_y_continuous(limits = c(-2000, 1000), breaks = seq(-2000, 1000, by = 500))

# print in a file
  pdf(file = "Figures/PCA_CTSS_merged_rep_PGCs_oocyte_embryo.pdf", height = 6, width = 8)
    print(p)
  dev.off()