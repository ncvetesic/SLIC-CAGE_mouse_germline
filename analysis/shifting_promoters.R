#/*==========================================================================#*/
#' ## Identifying shifting promoters
#/*==========================================================================#*/

# Preparation for shifting promoters - consensus clusters have to be called
  aggregateTagClusters(mESC_all_PGC_merge, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)
  cumulativeCTSSdistribution(mESC_all_PGC_merge, clusters = "consensusClusters")


# all vs E14_mESC - shifting promoters
  samples <- sampleLabels(mESC_all_PGC_merge)

# exclude E14_mESC
  samples <- samples[-1]

# call shifts
  all_vs_E14_mESC_shift.l <- lapply(samples, function(x) {
    scoreShift(mESC_all_PGC_merge, groupX = x, groupY = "E14_mESC", testKS = TRUE, useTpmKS = TRUE)
    shifting.promoters <- getShiftingPromoters(mESC_all_PGC_merge, tpmThreshold = 3, scoreThreshold = -Inf, fdrThreshold = 0.01)
    return(shifting.promoters) })
  names(all_vs_E14_mESC_shift.l) <- samples

# convert shifting list to GRanges
  all_vs_E14_mESC_shift.grl <- lapply(all_vs_E14_mESC_shift.l, function(x) {
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
    seqinfo(gr) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(gr)]
    return(gr)
  })

# save switching tss GRanges object
  saveRDS(all_vs_E14_mESC_shift.grl, "merged_replicates/intermediate_data/all_vs_E14_mESC_shift_grl.RDS")

# create shifting GRanges centred on mESC (group Y) - add distance separately
  library(BSgenome.Mmusculus.UCSC.mm10)
  
  all_vs_E14_mESC_shift.grl <- readRDS( "merged_replicates/intermediate_data/all_vs_E14_mESC_shift_grl.RDS")
  
# keep only standard chromosomes
  all_vs_E14_mESC_shift.grl <- lapply(all_vs_E14_mESC_shift.grl, function(x) {
    gr <- keepStandardChromosomes(x, pruning.mode = "coarse")
    return(gr)
  })
  
# create domTSS object centred on oocyte domTSS (group Y)
  domTSS_all_vs_E14_mESC_shift.grl <- lapply(all_vs_E14_mESC_shift.grl, function(x) {
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
    seqinfo(gr) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(gr)]
    
    gr <-keepStandardChromosomes(gr, pruning.mode = "coarse")
    return(gr)
  })
  
  saveRDS(domTSS_all_vs_E14_mESC_shift.grl, "merged_replicates/intermediate_data/domTSS_all_vs_E14_mESC_shift_grl.RDS")
  
  # create dominant TSS centred GRanges - all vs oocyte shift centred on X group (mESC)
  domTSS_all_vs_E14_mESC_shift_X.grl <- lapply(domTSS_all_vs_E14_mESC_shift.grl, function(x) {
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
    seqinfo(gr) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(gr)]
    
    keepStandardChromosomes(gr, pruning.mode = "coarse")
    return(gr)
  })
  
  for(i in 1:length(domTSS_all_vs_E14_mESC_shift.grl)) {
    domTSS_all_vs_E14_mESC_shift.grl[[i]]$domTSS_dist <- domTSS_all_vs_E14_mESC_shift.grl[[i]]$groupX.pos - domTSS_all_vs_E14_mESC_shift.grl[[i]]$groupY.pos
  }
  
  for(i in 1:length(domTSS_all_vs_E14_mESC_shift.grl)) {
    domTSS_all_vs_E14_mESC_shift.grl[[i]]$domTSS_dist[(as.character(strand(domTSS_all_vs_E14_mESC_shift.grl[[i]])) == "-")] <- -domTSS_all_vs_E14_mESC_shift.grl[[i]]$domTSS_dist[(as.character(strand(domTSS_all_vs_E14_mESC_shift.grl[[i]])) == "-")]
    
    domTSS_all_vs_E14_mESC_shift_X.grl[[i]]$domTSS_dist <- domTSS_all_vs_E14_mESC_shift.grl[[i]]$domTSS_dist 
  }
  saveRDS(domTSS_all_vs_E14_mESC_shift.grl, "merged_replicates/intermediate_data/domTSS_all_vs_E14_mESC_shift_grl.RDS")
  saveRDS(domTSS_all_vs_E14_mESC_shift_X.grl, "merged_replicates/intermediate_data/ddomTSS_all_vs_E14_mESC_shift_X_grl.RDS")
