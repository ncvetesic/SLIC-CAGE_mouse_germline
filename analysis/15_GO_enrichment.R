#/*==========================================================================#*/
#' ## Extended Data Figure 4A,B, 12A,B
#/*==========================================================================#*/
# GO enrichment - example code using shifting promoters, any input genes can be used

library(clusterProfiler)
# annotate shifting promoters
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)

  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  PGC_oocyte_embryo_vs_mESC_shift.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/PGC_oocyte_embryo_vs_mESC_shift_grl.RDS")

# use peak anno from ChipSeeker to annotate switching tag clusters
  peakAnno.l <- lapply(GRangesList(PGC_oocyte_embryo_vs_mESC_shift.grl), function(x) annotatePeak(x, TxDb = txdb,  annoDb = NULL, sameStrand = TRUE, verbose = FALSE))

# convert to GRanges 
  PGC_oocyte_embryo_vs_mESC_shift_anno.grl <- lapply(peakAnno.l, as.GRanges)

# save annotated GRanges object 
  saveRDS(PGC_oocyte_embryo_vs_mESC_shift_anno.grl, "../all_reps_PN6_2cell/intermediate_data/PGC_oocyte_embryo_vs_mESC_shift_anno_grl.RDS")

# load switching annotated GRanges 
  PGC_oocyte_embryo_vs_mESC_shift_anno.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/PGC_oocyte_embryo_vs_mESC_shift_anno_grl.RDS")

# load annotated tag clusters
  tc_PGCs_oocyte_embryo_anno.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/tc_PGCs_oocyte_embryo_anno_grl.RDS")

# create gene universe - all genes mouse annotated
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

  geneUniverse <- genes(txdb)$gene_id
  samples <- names(PGC_oocyte_embryo_vs_mESC_shift_anno.grl)

  geneSample.l <- lapply(PGC_oocyte_embryo_vs_mESC_shift_anno.grl, function(x) unique((x$geneId)))
  names(PGC_oocyte_embryo_vs_mESC_shift_anno.grl) <- samples


  ego_MF.l <- lapply(geneSample.l, function(x) 
    enrichGO(gene = x,
             universe      = geneUniverse,
             OrgDb         = org.Mm.eg.db,
             ont           = "MF",
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.05,
             readable      = FALSE,
             keyType = "ENTREZID"))

  ego_BP.l <- lapply(geneSample.l, function(x) 
    enrichGO(gene = x,
             universe      = geneUniverse,
             OrgDb         = org.Mm.eg.db,
             ont           = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.05,
             readable      = FALSE,
             keyType = "ENTREZID"))

  ego_CC.l <- lapply(geneSample.l, function(x) 
    enrichGO(gene = x,
             universe      = geneUniverse,
             OrgDb         = org.Mm.eg.db,
             ont           = "CC",
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.05,
             readable      = FALSE,
             keyType = "ENTREZID"))

# print output in pdf
  for (i in names(ego_BP.l)) {
    
    p <-  dotplot(ego_MF.l[[i]])
    pdf(paste0("merged_replicates/results/shifting_promoters/mESC/GO_enrichment/universe_mouse_annotated_all/MF/GO_MF_", i, ".pdf"), height = 4, width = 12)
    print(p)
    dev.off()
    
    p <-  dotplot(ego_BP.l[[i]])
    pdf(paste0("merged_replicates/results/shifting_promoters/mESC/GO_enrichment/universe_mouse_annotated_all/BP/GO_BP_", i, ".pdf"), height = 4, width = 12)
    print(p)
    dev.off()
    
    
    p <-  dotplot(ego_CC.l[[i]])
    pdf(paste0("merged_replicates/results/shifting_promoters/mESC/GO_enrichment/universe_mouse_annotated_all/CC/GO_CC_", i, ".pdf"), height = 4, width = 12)
    print(p)
    dev.off()
  }
