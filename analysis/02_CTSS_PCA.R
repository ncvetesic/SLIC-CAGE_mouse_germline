#/*==========================================================================#*/
#' ## Figure 1D
#/*==========================================================================#*/
# PCA plot at CTSS level using DeSEQ2

# need the raw TPM values, so I need to set the normalize method to "none", recreate consensus clusters and extract TPM
# copy myCAGEset
  rawCAGEset_merged_reps <- CAGEset_PGC_embryo_merged

# set raw TPM values
  normalizeTagCount(rawCAGEset_merged_reps, method = "none")

# clustering of CTSS (raw tag number) - beware everything with one read is clustered together
  clusterCTSS(object = rawCAGEset_merged_reps,
              threshold = 1,
              thresholdIsTpm = FALSE,
              nrPassThreshold = 1,
              method = "distclu",
              maxDist = 20,
              removeSingletons = TRUE, 
              keepSingletonsAbove = 5)

# cumulative distribution and quantile positions (raw tag number)
  cumulativeCTSSdistribution(rawCAGEset_merged_reps, clusters = "tagClusters")
  quantilePositions(rawCAGEset_merged_reps, 
                    clusters = "tagClusters", 
                    qLow = 0.1, 
                    qUp = 0.9)

# aggregate the clusters across the samples (raw tag number): remember these are raw reads in Tpm slot
  aggregateTagClusters(rawCAGEset_merged_reps, 
                       tpmThreshold = 5, 
                       qLow = 0.1, 
                       qUp = 0.9, 
                       maxDist = 100)

# save the CAGEset object (remember not normalized, only for Deseq2)
  saveRDS(rawCAGEset_merged_reps, file = "../all_reps_PN6_2cell/intermediate_data/rawCAGEset_merged_reps.RDS")

# load raw CAGE set
  rawCAGEset_merged_reps <- readRDS(file = "../all_reps_PN6_2cell/intermediate_data/rawCAGEset_merged_reps.RDS")

## Prepare the raw count table
# extract the raw counts
  raw.counts.df <- data.frame(consensusClustersTpm(rawCAGEset_merged_reps))

# extract information about consensus clusters (consensus coordinates)
  consensus.info <- consensusClusters(rawCAGEset_merged_reps)

# check dimensions of both tables
  dim(raw.counts.df)
  dim(consensus.info) # dimensions - number of rows are the same, I can just merge/combine the tables

# create identifiers 
  consensus.info$cons_clust_id <- paste("cid_", 1:nrow(consensus.info), sep = "")
  rownames(raw.counts.df) <- consensus.info$cons_clust_id

# combine the TPM table with the consensus cluster information
  count.consensus.info <- cbind(consensus.info[, -1], raw.counts.df)

# write table
  write.table(count.consensus.info, "../all_reps_PN6_2cell/intermediate_data/CountTable_raw_consensus.txt",
              col.names = TRUE,
              row.names = FALSE,
              sep = "\t",
              quote = FALSE)

# prepare DeSeq2 object
  library(DESeq2)

# create a count matrix for DESeq2
  count.matrix <- as.matrix(raw.counts.df)


#------- all merged data 
# info table/ design table / colData
# load raw CAGE set
  rawCAGEset_merged_reps <- readRDS(file = "../all_reps_PN6_2cell/intermediate_data/rawCAGEset_merged_reps.RDS")

## Prepare the raw count table
  raw.counts.df <- data.frame(consensusClustersTpm(rawCAGEset_merged_reps))

# prepare DeSeq2 object
  library(DESeq2)

# create a count matrix for DESeq2
  count.matrix <- as.matrix(raw.counts.df)



  samples <- colnames(count.matrix)
  info.df <- data.frame(condition = factor(x = c("mESC-E14", rep("PGC", times = 11), rep("oocyte", times  = 3), rep("embryo", times = 2)), 
                                           levels = c("mESC-E14", "PGC", "oocyte", "embryo")),
                        row.names = samples)

# create DESeqDataSet object - stores count matrix and design
  dds <- DESeqDataSetFromMatrix(countData = count.matrix,
                                colData = info.df,
                                design = ~ condition)

  col <- c("gold1", "darkslateblue", "dodgerblue")


# PCA preparation
  dds <- DESeq(dds)
  vsd <- vst(dds, blind = TRUE) # blind means don't take design into consideration

# save dds object PGC embryo oocyte  
  saveRDS(dds, "../all_reps_PN6_2cell/intermediate_data/dds_all_merg_no_TBP2.RDS")

# using top 10000 most variable tag clusters
  pcaData <- plotPCA(vsd, intgroup=c("condition"), ntop = 10000, returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

# plotting of PCA
  library(ggrepel)

# change names of samples and levels
  samples <- factor(c("mESC-E14", "E9.5", "E10.5", "E11.5", "E12.5F", "E12.5M", "E13.5F", "E13.5M", "E14.5F", "E14.5M", "E16.5F", "E16.5M", "P6", "P14", "MII", "2-cell", "4-cell"), levels = c("mESC-E14", "E9.5", "E10.5", "E11.5", "E12.5F", "E12.5M", "E13.5F", "E13.5M", "E14.5F", "E14.5M", "E16.5F", "E16.5M", "P6", "P14", "MII", "2-cell", "4-cell"))
  pcaData$name <- samples

# setting plotting colour
  col <- col_scheme

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
    scale_x_continuous(limits = c(-300, 500), breaks = seq(-300, 500, by = 100)) +
    scale_y_continuous(limits = c(-200, 200), breaks = seq(-200, 200, by = 100))



  pdf(file = "Figures/PCA_merged_rep_mESC_PGCs_oocyte_embryo.pdf", height = 8, width = 6)
    print(p)
  dev.off()