#/*==========================================================================#*/
#' ## Extended Data Figure 8D-G
#/*==========================================================================#*/
# SOM promoter class GO enrichment

# read in libraries
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# load som annotated GRanges
  consensus.clusters.som_promOnly.anno.grl <- readRDS("intermediate_data/consensus_clusters_som_promOnly_anno_Damir.9classes_grl.RDS")

# load annotated consensus clusters
  consensus.clusters.anno.gr <- readRDS("intermediate_data/consensus.clusters.anno.gr.RDS")

# create gene universe - all genes mouse annotated
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

  geneUniverse <- genes(txdb)$gene_id
  samples <- names(consensus.clusters.som_promOnly.anno.grl)
  
  geneSample.l <- lapply(consensus.clusters.som_promOnly.anno.grl, function(x) unique((x$geneId)))
  names(consensus.clusters.som_promOnly.anno.grl) <- samples


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
    pdf(paste0("merged_replicates/results/SOM_promOnly/GO_enrichment/universe_mouse_annotated_all/MF/GO_MF_", i, ".pdf"), height = 4, width = 10)
    print(p)
    dev.off()
    
    p <-  dotplot(ego_BP.l[[i]])
    pdf(paste0("merged_replicates/results/SOM_promOnly/GO_enrichment/universe_mouse_annotated_all/BP/GO_BP_", i, ".pdf"), height = 4, width = 10)
    print(p)
    dev.off()
    
    
    p <-  dotplot(ego_CC.l[[i]])
    pdf(paste0("merged_replicates/results/SOM_promOnly//GO_enrichment/universe_mouse_annotated_all/CC/GO_CC_", i, ".pdf"), height = 4, width = 10)
    print(p)
    dev.off()
  }

# save results in txt files
  ego_CC_txt.l <- list()
  ego_BP_txt.l <- list()
  ego_MF_txt.l <- list()

# CC mESC add info on SOM class
  for (i in 1:length(ego_CC.l)) {
    ego_CC_txt.l[[i]] <- as.data.frame(ego_CC.l[[i]])
    ego_CC_txt.l[[i]]$SOM = rep(names(ego_CC.l[i]), times = nrow(ego_CC.l[[i]]))
  }

  write.csv(do.call(rbind, ego_CC_txt.l), file = "merged_replicates/results/SOM_promOnly/GO_enrichment/universe_mouse_annotated_all/CC/GO_CC.csv")

# BP add info on SOM class
  for (i in 1:length(ego_BP.l)) {
    ego_BP_txt.l[[i]] <- as.data.frame(ego_BP.l[[i]])
    ego_BP_txt.l[[i]]$SOM = rep(names(ego_BP.l[i]), times = nrow(ego_BP.l[[i]]))
  }

  write.csv(do.call(rbind, ego_BP_txt.l), file = "merged_replicates/results/SOM_promOnly/GO_enrichment/universe_mouse_annotated_all/BP/GO_BP.csv")    

# MF add info on SOM class
  for (i in 1:length(ego_MF.l)) {
    ego_MF_txt.l[[i]] <- as.data.frame(ego_MF.l[[i]])
    ego_MF_txt.l[[i]]$SOM = rep(names(ego_MF.l[i]), times = nrow(ego_MF.l[[i]]))
  }

  write.csv(do.call(rbind, ego_MF_txt.l), file = "merged_replicates/results/SOM_promOnly/GO_enrichment/universe_mouse_annotated_all/MF/GO_MF_csv")   