# note this is example processing - not all data is loaded at this stage, but all parameters are the same
# in power-law normalisation use alpha as suggested from power-law distribution

# Install old CAGEr - 1.20.0 version was used

  devtools::install_version("CAGEr", version = "1.20.0", repos = "http://bioconductor.org/packages/3.6/bioc/")
  library(CAGEr)


# check CAGEr version - has to be 1.16.0 or < 1.20.0
  packageVersion("CAGEr")

  library(BSgenome.Mmusculus.UCSC.mm10)

# set input directory where mapped bam files E14 nantiCAGE and get paths to sorted.bam files:
  inputDir <- ("../../../Hajkova_collaboration/J1_TKO_comparison/data/170609_D00467_0254_Aca386anxx_J1_TKO/")
  paths <- list.files(inputDir, full.names = TRUE)
  pathsToInputFiles0<- paths[grep(paths, pattern = "E14_sample8_hiseq.sorted.bam$")]

  inputDir <- "../../mapped/merged_lanes/"
  paths <- list.files(inputDir, full.names = TRUE)
  pathsToInputFiles1 <- paths[grep(paths, pattern = ".bam$")]
  pathsToInputFiles1 <- pathsToInputFiles1[c(18, 17, 3:1, 4, 7, 5, 8, 6, 11, 9, 12, 10, 13:16, 19:21)]
  
  inputDir <- "../../mapped/rep3_2019_April/"
  paths <- list.files(inputDir, full.names = TRUE)
  pathsToInputFiles2 <- paths[grep(paths, pattern = ".sorted.bam$")]
  pathsToInputFiles2 <- pathsToInputFiles2[c(19, 18, 15, 16, 17)]


# combine all paths to inputFiles
  pathsToInputFiles <- c(pathsToInputFiles0, pathsToInputFiles1, pathsToInputFiles2)

# create myCAGEset object
  mESC_all_PGC_CAGE <- new("CAGEset", genomeName = "BSgenome.Mmusculus.UCSC.mm10",
                           inputFilesType = "bam",
                           inputFiles = pathsToInputFiles, 
                           sampleLabels = c("E14_mESC", "E9_5_r1", "E9_5_r2", 
                                            "E10_5_r1", "E10_5_r2", "E10_5_r3",
                                            "E11_5_r1",
                                            "E12_5F_r1", "E12_5F_r2",
                                            "E12_5M_r1", "E12_5M_r2",
                                            "E13_5F_r1", "E13_5F_r2",
                                            "E13_5M_r1", "E13_5M_r2",
                                            "E16_5F_r1", "E16_5F_r2",
                                            "E16_5M_r1", "E16_5M_r2",
                                            "oocyte_r1", "oocyte_r2", "oocyte_r3",
                                            "PN7_r1", "PN7_TBP2_KO_r1",
                                            "PN14_r1", "PN14_r2",
                                            "PN14_TBP2_KO_r1"))


  getCTSS(mESC_all_PGC_CAGE)

# save CAGE unmerged object
  saveRDS(mESC_all_PGC_CAGE, "intermediate_data/mESC_all_PGC_CAGE.RDS")


# Correlation of the raw data
# plot correlation of raw data
  setwd("qc")
  corr.m <- plotCorrelation(mESC_all_PGC_CAGE, samples = "all", method = "pearson",
                            what = "CTSS", values = "raw", tagCountThreshold = 1, applyThresholdBoth = TRUE)

# Normalization of the libraries
  #check the sizes of libraries
    librarySizes(mESC_all_PGC_CAGE) 

    setwd("qc/")
    plotReverseCumulatives(mESC_all_PGC_CAGE, fitInRange = c(10, 10000), onePlot = TRUE, values = "raw")

#normalize using power-law distribution, alpha = 1.14
  normalizeTagCount(mESC_all_PGC_CAGE, method = "powerLaw", fitInRange = c(10, 10000), alpha = 1.14, T = 10^6)

#plot normalized reverse cumulative number of CAGE tags per CTSS
  plotReverseCumulatives(mESC_all_PGC_CAGE, fitInRange = c(10, 10000), onePlot = TRUE, values = "normalized")

# Correlation of the normalized data
  # plot correlation of raw data
    setwd("qc")
    corr.m <- plotCorrelation(mESC_all_PGC_CAGE, samples = "all", method = "pearson",
                              what = "CTSS", values = "normalized", tagCountThreshold = 1, applyThresholdBoth = TRUE)

# Create tracks - export CAGE signal (CTSS) to bedgraph
    setwd("tracks/")
    exportCTSStoBedGraph(mESC_all_PGC_CAGE, values = "normalized", format = "bedGraph", oneFile = FALSE)

    
# Cluster CTSSs into tag clusters
    clusterCTSS(object = mESC_all_PGC_CAGE, threshold = 1, thresholdIsTpm = TRUE,
                nrPassThreshold = 1, method = "distclu", maxDist = 20,
                removeSingletons = TRUE, keepSingletonsAbove = 5)
    
# Promoter width
    setwd("results")
    cumulativeCTSSdistribution(mESC_all_PGC_CAGE, clusters = "tagClusters")
    quantilePositions(mESC_all_PGC_CAGE, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
    plotInterquantileWidth(mESC_all_PGC_CAGE, clusters = "tagClusters", tpmThreshold = 3, qLow = 0.1, qUp = 0.9)
    
# Extract tag clusters - from non-merged replicates
    samples <- sampleLabels(mESC_all_PGC_CAGE)
    quantilePositions(mESC_all_PGC_CAGE, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
    
    tc_PGC_reps.l <- list()
    
    for (i in 1:length(samples)) {
      tc_PGC_reps.l[[i]] <- tagClusters(mESC_all_PGC_CAGE, sample = samples[i], returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
      names(tc_PGC_reps.l)[[i]] <- samples[i]
    }
    
    #check the number of tag clusters/ promoters
    sapply(tc_PGC_reps.l, nrow)
    
    #save tag clusters for later use
    for (i in 1:length(samples)) {  
      write.table(tc_PGC_reps.l[[i]], 
                  paste("intermediate_data/tc_PGC_reps_", samples[i], ".txt", sep = ""), 
                  col.names = TRUE,
                  row.names = FALSE,
                  quote = FALSE,
                  sep = "\t")
    }

    
# Convert tag cluster dataframe into GRanges
    tc_PGC_reps.grl <- list()
    
    for (i in 1:length(samples)) {
      tc_PGC_reps.grl[[i]] <- GRanges(seqnames = tc_PGC_reps.l[[i]]$chr,
                                      ranges = IRanges(start = tc_PGC_reps.l[[i]]$start, 
                                                       end = tc_PGC_reps.l[[i]]$end),
                                      strand = tc_PGC_reps.l[[i]]$strand,
                                      nr_ctss = tc_PGC_reps.l[[i]]$nr_ctss,
                                      dominant_ctss = tc_PGC_reps.l[[i]]$dominant_ctss,
                                      tpm = tc_PGC_reps.l[[i]]$tpm,
                                      tpm.dominant_ctss = tc_PGC_reps.l[[i]]$tpm.dominant_ctss,
                                      q_0.1 = tc_PGC_reps.l[[i]]$q_0.1,
                                      q_0.9 = tc_PGC_reps.l[[i]]$q_0.9,
                                      interquantile_width = tc_PGC_reps.l[[i]]$interquantile_width)
      
      seqlevels(tc_PGC_reps.grl[[i]]) <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
      
  # attach chromosome sequence lengths
    seqlengths(tc_PGC_reps.grl[[i]]) <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
        
  # attach genome information
    genome(tc_PGC_reps.grl[[i]]) <- genome(BSgenome.Mmusculus.UCSC.mm10)
      
  # drop the scaffolds - keep standard chromosomes
    tc_PGC_reps.grl[[i]] <- keepStandardChromosomes(tc_PGC_reps.grl[[i]], pruning.mode = "coarse", species = "Mus_musculus")
    }
    
# attach sample names
  names(tc_PGC_reps.grl) <- samples
    
# save tag cluster granges object for later
  saveRDS(tc_PGC_reps.grl, file ="intermediate_data/tc_PGC_reps_grl.RDS")

# Create GRanges object from dominant TSS positions (merged datasets per lane)
  domTSS_PGC_reps.grl <- list()
    
    for (i in 1:length(samples)) {
      domTSS_PGC_reps.grl[[i]] <- GRanges(seqnames = tc_PGC_reps.l[[i]]$chr,
                                          ranges = IRanges(start = tc_PGC_reps.l[[i]]$dominant_ctss, 
                                                           end = tc_PGC_reps.l[[i]]$dominant_ctss),
                                          strand = tc_PGC_reps.l[[i]]$strand,
                                          nr_ctss = tc_PGC_reps.l[[i]]$nr_ctss,
                                          dominant_ctss = tc_PGC_reps.l[[i]]$dominant_ctss,
                                          tpm = tc_PGC_reps.l[[i]]$tpm, 
                                          tpm.dominant_ctss = tc_PGC_reps.l[[i]]$tpm.dominant_ctss, 
                                          q_0.1 = tc_PGC_reps.l[[i]]$q_0.1,
                                          q_0.9 = tc_PGC_reps.l[[i]]$q_0.9,
                                          interquantile_width = tc_PGC_reps.l[[i]]$interquantile_width)
      
      #attach genome information
      seqinfo(domTSS_PGC_reps.grl[[i]]) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(domTSS_PGC_reps.grl[[i]])]
    }
    
# attach sample names
  names(domTSS_PGC_reps.grl) <- samples
    
# save tag cluster granges object for later
  saveRDS(domTSS_PGC_reps.grl, file ="intermediate_data/domTSS_PGC_reps_grl.RDS")

