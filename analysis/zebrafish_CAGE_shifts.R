#/*==========================================================================#*/
#' ## Zebrafish CAGE data
#/*==========================================================================#*/


# import Zebrafish developmental CAGE
  install.packages("../../extra_data/ZebrafishDevelopmentalCAGE_0.99.0.tar.gz", repos = NULL, type="source")
  library(ZebrafishDevelopmentalCAGE)
  library(BSgenome.Drerio.UCSC.danRer7)
  data(ZebrafishSamples)

# import zebrafish samples = egg - maternal and prim6 - somatic, after ZGA
  zebrafishCAGEset <- importPublicData(source = "ZebrafishDevelopment", dataset = "ZebrafishCAGE", 
                                       group = "development", sample = c("zf_unfertilized_egg", "zf_prim6"))

# plot correlation of raw data
  setwd("zebrafish/")
  corr.m <- plotCorrelation(zebrafishCAGEset, samples = "all", method = "pearson",
                            what = "CTSS", values = "raw", tagCountThreshold = 1, applyThresholdBoth = TRUE)


# normalization of the libraries
## check the sizes of libraries
  librarySizes(zebrafishCAGEset) 

  plotReverseCumulatives(zebrafishCAGEset, fitInRange = c(5, 10000), onePlot = TRUE, values = "raw")

# normalize using power-law distribution, alpha = 1.24
  normalizeTagCount(zebrafishCAGEset, method = "powerLaw", fitInRange = c(5, 10000), alpha = 1.24, T = 10^6)

# plot normalized reverse cumulative number of CAGE tags per CTSS
  plotReverseCumulatives(zebrafishCAGEset, fitInRange = c(5, 10000), onePlot = TRUE, values = "normalized")

## correlation of the normalized data
# plot correlation of raw data
corr.m <- plotCorrelation(zebrafishCAGEset, samples = "all", method = "pearson",
                          what = "CTSS", values = "normalized", tagCountThreshold = 1, applyThresholdBoth = TRUE)

# cluster CTSSs into tag clusters
  clusterCTSS(object = zebrafishCAGEset, threshold = 1, thresholdIsTpm = TRUE,
              nrPassThreshold = 1, method = "distclu", maxDist = 20,
              removeSingletons = TRUE, keepSingletonsAbove = 5)

# export normalised bedgraph tracks
  setwd("zebrafish/")
  exportCTSStoBedGraph(zebrafishCAGEset, values = "normalized", format = "bedGraph", oneFile = FALSE)

# promoter width
  setwd("zebrafish")
  cumulativeCTSSdistribution(zebrafishCAGEset, clusters = "tagClusters")
  quantilePositions(zebrafishCAGEset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
  plotInterquantileWidth(zebrafishCAGEset, clusters = "tagClusters", tpmThreshold = 3, qLow = 0.1, qUp = 0.9)

# save zebrafish CAGE normalised object
  saveRDS(zebrafishCAGEset, "zebrafishCAGEset.RDS")

# extract tag clusters
  samples <- sampleLabels(zebrafishCAGEset)

  tc_fish.l <- list()
  
  for (i in 1:length(samples)) {
    tc_fish.l[[i]] <- tagClusters(zebrafishCAGEset, sample = samples[i], returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
    names(tc_fish.l)[[i]] <- samples[i]
  }

# check the number of tag clusters/ promoters
  sapply(tc_fish.l, nrow)

# save tag clusters for later use
  for (i in 1:length(samples)) {  
    write.table(tc_fish.l[[i]], 
                paste("zebrafish/tc_fish_", samples[i], ".txt", sep = ""), 
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
  }

# convert tag cluster dataframe into GRanges
  tc_fish.grl <- list()

  for (i in 1:length(samples)) {
    tc_fish.grl[[i]] <- GRanges(seqnames = tc_fish.l[[i]]$chr,
                                ranges = IRanges(start = tc_fish.l[[i]]$start, 
                                                 end = tc_fish.l[[i]]$end),
                                strand = tc_fish.l[[i]]$strand,
                                nr_ctss = tc_fish.l[[i]]$nr_ctss,
                                dominant_ctss = tc_fish.l[[i]]$dominant_ctss,
                                tpm = tc_fish.l[[i]]$tpm,
                                tpm.dominant_ctss = tc_fish.l[[i]]$tpm.dominant_ctss,
                                q_0.1 = tc_fish.l[[i]]$q_0.1,
                                q_0.9 = tc_fish.l[[i]]$q_0.9,
                                interquantile_width = tc_fish.l[[i]]$interquantile_width)
  
  # attach genome information
  seqinfo(tc_fish.grl[[i]]) <- seqinfo(BSgenome.Drerio.UCSC.danRer7)[seqlevels(tc_fish.grl[[i]])]
  
  # drop the scaffolds - keep standard chromosomes
    tc_fish.grl[[i]] <- keepStandardChromosomes(tc_fish.grl[[i]], pruning.mode = "coarse")
  }

# attach sample names
  names(tc_fish.grl) <- samples

# save tag cluster granges object for later
  saveRDS(tc_fish.grl, file ="zebrafish/tc_fish_grl.RDS")

# prepare for shifting promoters
  aggregateTagClusters(zebrafishCAGEset, tpmThreshold = 3, qLow = 0.1, qUp = 0.9, maxDist = 100)
  cumulativeCTSSdistribution(zebrafishCAGEset, clusters = "consensusClusters")

# save cage set with consensus called
  saveRDS(zebrafishCAGEset, "zebrafish/zebrafishCAGEset.RDS")

  samples <- sampleLabels(zebrafishCAGEset)[1]

# shifting promoters between zf_egg and zf_prim6
  zf_egg_vs_prim6_shift.l <- lapply(samples, function(x) {
    scoreShift(zebrafishCAGEset, groupX = x, groupY = "zf_prim6", testKS = TRUE, useTpmKS = TRUE)
    shifting.promoters <- getShiftingPromoters(zebrafishCAGEset, tpmThreshold = 3, scoreThreshold = -Inf, fdrThreshold = 0.01)
    return(shifting.promoters)
    })
  
  names(zf_egg_vs_prim6_shift.l) <- samples

# convert to GRanges
  zf_egg_vs_prim6_shift.grl <- lapply(zf_egg_vs_prim6_shift.l, function(x) {
    gr <- GRanges(seqnames = x$chr, 
                  ranges = IRanges(start = x$start,
                                   end = x$end),
                  strand = x$strand,
                  shifting.score = x$shifting.score,
                  groupX.pos = x$groupX.pos,
                  groupY.pos = x$groupY.pos,
                  groupX.tpm = x$groupX.tpm,
                  groupY.tpm = x$groupY.tpm,
                  pvalue.KS = x$pvalue.KS,
                  fdr.KS = x$fdr.KS)
    seqinfo(gr) <- seqinfo(BSgenome.Drerio.UCSC.danRer7)[seqlevels(gr)]
    return(gr)
  })

# save switching tss GRanges object
  saveRDS(zf_egg_vs_prim6_shift.grl, "zebrafish/zf_egg_vs_prim6_shift_grl.RDS")

# create GRanges centred on X or Y group dominant TSS
  zf_egg_vs_prim6_shift.grl <- readRDS("zebrafish/zf_egg_vs_prim6_shift_grl.RDS")

# create domTSS object centred on oocyte domTSS (group x)
  domTSS_egg_p6_shift_fish_egg_centred.grl <- lapply(zf_egg_vs_prim6_shift.grl, function(x) {
    gr <- GRanges(seqnames = seqnames(x), 
                  ranges = IRanges(start = x$groupX.pos,
                                   end = x$groupX.pos),
                  strand = strand(x),
                  start_tc = start(x),
                  end_tc = end(x),
                  groupX.pos = x$groupX.pos,
                  groupY.pos = x$groupY.pos,
                  groupX.tpm = x$groupX.tpm,
                  groupY.tpm = x$groupY.tpm,
                  shifting.score = x$shifting.score,
                  pvalue.KS = x$pvalue.KS,
                  fdr.KS = x$fdr.KS)
    seqinfo(gr) <- seqinfo(BSgenome.Drerio.UCSC.danRer7)[seqlevels(gr)]
  
    gr <-keepStandardChromosomes(gr, pruning.mode = "coarse")
    return(gr)
  })


# create domTSS object centred on prim6 domTSS (group y)
domTSS_egg_p6_shift_fish_prim6_centred.grl <- lapply(zf_egg_vs_prim6_shift.grl, function(x) {
  gr <- GRanges(seqnames = seqnames(x), 
                ranges = IRanges(start = x$groupY.pos,
                                 end = x$groupY.pos),
                strand = strand(x),
                start_tc = start(x),
                end_tc = end(x),
                groupX.pos = x$groupX.pos,
                groupY.pos = x$groupY.pos,
                groupX.tpm = x$groupX.tpm,
                groupY.tpm = x$groupY.tpm,
                shifting.score = x$shifting.score,
                pvalue.KS = x$pvalue.KS,
                fdr.KS = x$fdr.KS)
  seqinfo(gr) <- seqinfo(BSgenome.Drerio.UCSC.danRer7)[seqlevels(gr)]
  
  keepStandardChromosomes(gr, pruning.mode = "coarse")
  return(gr)
})

# attach distances between dominant TSSs  
  for(i in 1:length(domTSS_egg_p6_shift_fish_egg_centred.grl)) {
    domTSS_egg_p6_shift_fish_egg_centred.grl[[i]]$domTSS_dist <- domTSS_egg_p6_shift_fish_egg_centred.grl[[i]]$groupX.pos - domTSS_egg_p6_shift_fish_egg_centred.grl[[i]]$groupY.pos
  }

  for(i in 1:length(domTSS_egg_p6_shift_fish_egg_centred.grl)) {
    domTSS_egg_p6_shift_fish_egg_centred.grl[[i]]$domTSS_dist[(as.character(strand(domTSS_egg_p6_shift_fish_egg_centred.grl[[i]])) == "-")] <- -domTSS_egg_p6_shift_fish_egg_centred.grl[[i]]$domTSS_dist[(as.character(strand(domTSS_egg_p6_shift_fish_egg_centred.grl[[i]])) == "-")]
    domTSS_egg_p6_shift_fish_prim6_centred.grl[[i]]$domTSS_dist <- domTSS_egg_p6_shift_fish_egg_centred.grl[[i]]$domTSS_dist 
  }

# save domTSS centred zebrafish objects
  saveRDS(domTSS_egg_p6_shift_fish_egg_centred.grl, "zebrafish/domTSS_egg_p6_shift_fish_egg_centred_grl.RDS")
  saveRDS(domTSS_egg_p6_shift_fish_prim6_centred.grl, "zebrafish/domTSS_egg_p6_shift_fish_prim6_centred_grl.RDS")