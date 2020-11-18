#/*==========================================================================#*/
#' ## Extended Data Figure 15A
#/*==========================================================================#*/
# dominant CTSS classification using SOMS - Self Organising Maps

# load libraries
  library(kohonen)
  library(tibble)
  library(stringr)
  library(forcats)

# load annotated consensus clusters
  consensus.clusters.anno.gr <- readRDS("intermediate_data/consensus.clusters.anno.gr.RDS")
  consensus.clusters.anno.gr <- consensus.clusters.anno.gr[[1]] # for some reason its saved as a list with 1 element, so this is just "unlisting"

# select promoter regions only
  promoters <- c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")
  consensus.clusters.anno.prom.gr <- consensus.clusters.anno.gr[consensus.clusters.anno.gr$annotation %in% promoters]


# CTSS sample selection - exclude bad samples and thse that just might interfere -  E9.5, P7 TBP2 KO, and embryo, and focus only on female as there is transition to oocyte too
  samples_sel <- c("E14_mESC", "E10_5", "E11_5", "E12_5F", "E13_5F", 
                   "E14_5F", "E16_5F", "PN6", "PN14", "oocyte")  
  samples <- samples_sel

#--to extract dominant TSSs i need a dominant TSS object - sample specific in a list, concatenate to one .gr object and overlap with CTSS.gr object
# extract normalized CTSSs
  CTSSnorm <- CTSSnormalizedTpm(CAGEset_PGC_embryo_merged)

# add CTSS number - so I can refer to them later to classify into SOM classes
  CTSSnorm$CTSSid <- 1:nrow(CTSSnorm)

# create GRanges to overlap with consensus cluster promoters
  CTSSnorm.gr <- GRanges(seqnames = CTSSnorm$chr,
                         ranges = IRanges(start = CTSSnorm$pos,
                                          end = CTSSnorm$pos),
                         strand = CTSSnorm$strand)
  mcols(CTSSnorm.gr) <- CTSSnorm[, samples_sel]
  CTSSnorm.gr$CTSSid <- CTSSnorm$CTSSid
  seqinfo(CTSSnorm.gr) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(CTSSnorm.gr)]

# overlap CTSS norm granges with consensus cluster granges to select only promoter level ctss
  CTSSnorm_prom.gr <- subsetByOverlaps(CTSSnorm.gr, consensus.clusters.anno.prom.gr)

# no before
  length(CTSSnorm.gr)

# no after filtering
  length(CTSSnorm_prom.gr)

# import domTSS.grl object - sample specific
  domTSS_PGCs_oocyte_embryo_anno.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/domTSS_PGCs_oocyte_embryo_anno_grl.RDS")

# select samples
  domTSS_PGCs_oocyte_embryo_anno_sel.grl <- domTSS_PGCs_oocyte_embryo_anno.grl[samples_sel]

# add sample mcols
  domTSS_PGCs_oocyte_embryo_anno_sel.grl <- lapply(1:length(domTSS_PGCs_oocyte_embryo_anno_sel.grl), function(x) {
    domTSS_PGCs_oocyte_embryo_anno_sel.grl[[x]]$sample <- samples_sel[x]
    return(domTSS_PGCs_oocyte_embryo_anno_sel.grl[[x]])
  })

# convert domTSS object to a single gr - beware unlisting will not work when you have same same elements
  domTSS_PGCs_oocyte_embryo_anno_sel.gr <- do.call(c, unlist(domTSS_PGCs_oocyte_embryo_anno_sel.grl, recursive = T, use.names = T))

# overlap with CTSS.gr object
  domCTSSnorm_prom.gr <- subsetByOverlaps(CTSSnorm_prom.gr, domTSS_PGCs_oocyte_embryo_anno_sel.gr)

# check that nothing is duplicated
  length(unique(domCTSSnorm_prom.gr$CTSSid)) # same length as length of .gr object

# convert to dataframe
  domCTSSnorm_prom.mat <- as.matrix(mcols(domCTSSnorm_prom.gr))

# add CTSSid as rownames 
  rownames(domCTSSnorm_prom.mat) <- domCTSSnorm_prom.mat[, "CTSSid"]

# omit the CTSSid column - last column
  domCTSSnorm_prom.mat <- domCTSSnorm_prom.mat[, -ncol(domCTSSnorm_prom.mat)]

## CAGEr normalization 1 tpm CTSS in at least 2 samples - to make it more robust
  norm.for.som <- function(tpm.mx, tpmThreshold = 1, nrPassThreshold = 2){
    sample.labels <- colnames(tpm.mx)
    nr.pass.threshold <- apply(tpm.mx, 1, function(x) {sum(x >= tpmThreshold)})
    idx <- nr.pass.threshold >= min(nrPassThreshold, length(sample.labels))
    ## scale divides the column by the scaling factor
    ## the scaling factor is calculate by applying sqrt(sum(x$ESCday0^2)/(n-1)), where n = number of rows with non-missing values
    ## non-negative log(tpm) are scaled by 
    ## by transposing the matrix with t, the scaling factor is calculated for each element over all the samples
    matForSom <- t(base::scale(t(log(tpm.mx+1)), center=F)) 
    matForSom <- matForSom[idx,]
    return(matForSom)
  }

  motForSom_Cons <- norm.for.som(domCTSSnorm_prom.mat)
  motForSom_Cons <- na.omit(motForSom_Cons)

# after all filtering I have  46571 CTSSs across 10 samples (and had 79888)
  dim(motForSom_Cons)

## som here is needed to cluster each element based on normalised score - set seed to have it reproducible
  set.seed(7)
  openSom_cons <- kohonen::som(motForSom_Cons,
                               grid = somgrid(5, 5, 'hexagonal'),
                               mode = 'pbatch',
                               cores = 4)
  
  somDf_cons <- as.data.frame(motForSom_Cons)
  somDf_cons <- rownames_to_column(somDf_cons, var = "elID")
  somDf_cons$elID <- as.numeric(somDf_cons$elID)
  somDf_cons$som <- openSom_cons$unit.classif

# windsorize signal - wont windsorize at this point
#somDf_cons_wins <- somDf_cons %>% dplyr::group_by(som) %>% dplyr::select(-c(som, elID)) %>% purrr::map(function(x) Winsorize(x)) %>% as.data.frame() 
  somDf_cons_wins <- somDf_cons
  somDf_cons_wins$elID <- somDf_cons$elID

## for consensusClusteres
  library(tidyr)
  somDf_cons <- gather(somDf_cons, "sample", "value", -c("elID", "som")) # transposes the columns in rows
  somDf_cons_wins <- gather(somDf_cons_wins, "sample", "value", -c("elID", "som")) # transposes the columns in rows


  somStats_cons <- somDf_cons %>% group_by(som, sample) %>% dplyr::summarize(avg = median(value))
  somStats_cons_wins <- somDf_cons_wins %>% group_by(som, sample) %>% dplyr::summarize(avg = median(value))
  
  somStats_cons$sample <- factor(somStats_cons$sample, levels = samples)
  somStats_cons_wins$sample <- factor(somStats_cons_wins$sample, levels = samples)
  
  somDf_cons$sample <- factor(somDf_cons$sample, levels = samples)
  somDf_cons_wins$sample <- factor(somDf_cons_wins$sample, levels = samples)

# count the number of clusters per SOM
  som_count <- somDf_cons %>% dplyr::group_by(som, sample) %>% dplyr::select(elID) %>% dplyr::tally()
  som_count_summary <- data.frame("som "= unique(som_count$som), "number" = unique(som_count$n))
  
  som_count_wins <- somDf_cons_wins %>% dplyr::group_by(som, sample) %>% dplyr::select(elID) %>% dplyr::tally()
  som_count_summary_wins <- data.frame("som "= unique(som_count_wins$som), "number" = unique(som_count_wins$n))


# add number to som class
  somDf_cons$som_num <- paste0(somDf_cons$som_num, " (", som_count_summary[somDf_cons$som, ]$number, ")")
  somStats_cons$som_num <- paste0(somStats_cons$som, " (", som_count_summary[somStats_cons$som, ]$number, ")")
  
  somDf_cons_wins$som_num <- paste0(somDf_cons_wins$som, " (", som_count_summary_wins[somDf_cons_wins$som, ]$number, ")")
  somStats_cons_wins$som_num <- paste0(somStats_cons_wins$som, " (", som_count_summary_wins[somStats_cons_wins$som, ]$number, ")")
  
  col <- c(col_scheme[-c(2, 6, 8, 10, 12, 14, 17, 18)])
  names(col) <- samples

# change levels for som_num so panels are printed in the right order
  somDf_cons_wins$som_num <- factor(somDf_cons_wins$som_num, 
                                    levels = c("1 (1429)", "2 (747)", "3 (1121)", "4 (1859)", "5 (3769)",
                                               "6 (1571)", "7 (2012)", "8 (2098)", "9 (2026)", "10 (1978)",
                                               "11 (1151)", "12 (4038)", "13 (3446)", "14 (1311)", "15 (1547)",
                                               "16 (1083)", "17 (2582)", "18 (2766)", "19 (1918)", "20 (1198)",
                                               "21 (1312)", "22 (1116)", "23 (1128)", "24 (1675)", "25 (1690)"))
  somStats_cons_wins$som_num <- factor(somStats_cons_wins$som_num,
                                       levels = c("1 (1429)", "2 (747)", "3 (1121)", "4 (1859)", "5 (3769)",
                                                  "6 (1571)", "7 (2012)", "8 (2098)", "9 (2026)", "10 (1978)",
                                                  "11 (1151)", "12 (4038)", "13 (3446)", "14 (1311)", "15 (1547)",
                                                  "16 (1083)", "17 (2582)", "18 (2766)", "19 (1918)", "20 (1198)",
                                                  "21 (1312)", "22 (1116)", "23 (1128)", "24 (1675)", "25 (1690)"))

##plot
  library(ggplot2)
  p <-  ggplot(somDf_cons_wins, aes(x = sample, y = value, fill = sample), inherit.aes = FALSE) + 
    geom_violin(alpha = .8, trim = TRUE, lwd = 0.25) +
    scale_fill_manual(values = col) +
    geom_point(data = somStats_cons_wins, aes(x = sample, y = avg), inherit.aes = FALSE, size = 0.5) + 
    geom_line(data = somStats_cons_wins, aes(group = 1, x = sample, y = avg), inherit.aes = FALSE, lwd = 0.25) +
    scale_color_manual(values = col_scheme) +
    theme_bw() +
    theme(text = element_text(size = 16),
          axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 1, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          strip.background =element_rect(fill="white")) +
    coord_cartesian(ylim=(c(0, 4)))
  
  pdf("Figures/SOM_domCTSSs_promOnly_F_5_5.pdf", height = 12, width = 14)
    p + facet_wrap(~ som_num)
  dev.off() 


# ----- extract SOM genes ----- #
# consensus cluster ID corresponds to elID
## extract consensus clusters
# split consensus cluster IDs according to SOM class - so consensus IDs will be duplicated because I have multiple samples and consensus clusters are same in each if in som class
# soms
  som_df.l <- split(somDf_cons_wins, list(somDf_cons_wins$som))
  
  # extract consensus cluster IDs - all are the same regardless of sample, so I will extract 25 SOM ids from mESC
  consID_som.l <- list()
  for(i in 1:25) {
    consID_som.l[[i]] <- as.numeric(unique(som_df.l[[i]]$elID))
  }
  names <- 1:25

# divide consensus clusters according to ID
  domCTSSs_SOM25_PGC_F_oocyte_promOnly.l <- list()

  for(i in 1:length(consID_som.l)) {
    domCTSSs_SOM25_PGC_F_oocyte_promOnly.l[[i]] <- CTSSnorm.gr[CTSSnorm.gr$CTSSid %in% consID_som.l[[i]], ]
  }

# convert go GRangesList
  domCTSSs_SOM25_PGC_F_oocyte_promOnly.grl <- GRangesList(domCTSSs_SOM25_PGC_F_oocyte_promOnly.l)
  rm(domCTSSs_SOM25_PGC_F_oocyte_promOnly.l)
  names(domCTSSs_SOM25_PGC_F_oocyte_promOnly.grl) <- paste0("som", 1:25)

# save SOM classes
  saveRDS(domCTSSs_SOM25_PGC_F_oocyte_promOnly.grl, "intermediate_data/domCTSSs_SOM25_PGC_F_oocyte_promOnly.grl.RDS")