#/*==========================================================================#*/
#' ## Extended Data Figure 14E, F
#/*==========================================================================#*/
# heatmap patterns in E16.5F vs mESC shifting single and double shift promoters


# load libraries
  library(heatmaps)
  library(BSgenome.Mmusculus.UCSC.mm10)

# objects needed
  double_shift_E16_5F_vs_mESC_housekeep_maternal.gr
  double_shift_E16_5F_vs_mESC_housekeep_somatic.gr
  single_shift_E16_5F_vs_mESC_housekeep.gr

# read nucleosome pwm object
  nucleosome_pwm <- readRDS("intermediate_data/nucleosome_pwm.RDS")

# --- calculate domTSS distances ---X is E16_5F and Y is mESC #
# double shift somatic
  double_shift_E16_5F_vs_mESC_housekeep_somatic.gr$domTSS_dist_E16_5FvsmESC <- double_shift_E16_5F_vs_mESC_housekeep_somatic.gr$groupX.pos - double_shift_E16_5F_vs_mESC_housekeep_somatic.gr$groupY.pos

  double_shift_E16_5F_vs_mESC_housekeep_somatic.gr$domTSS_dist_E16_5FvsmESC[as.character(strand(double_shift_E16_5F_vs_mESC_housekeep_somatic.gr)) == "-"] <- 
    -double_shift_E16_5F_vs_mESC_housekeep_somatic.gr$domTSS_dist_E16_5FvsmESC[as.character(strand(double_shift_E16_5F_vs_mESC_housekeep_somatic.gr)) == "-"]

# double shift maternal
  double_shift_E16_5F_vs_mESC_housekeep_maternal.gr$domTSS_dist_E16_5FvsmESC <- double_shift_E16_5F_vs_mESC_housekeep_maternal.gr$groupX.pos - double_shift_E16_5F_vs_mESC_housekeep_maternal.gr$groupY.pos

  double_shift_E16_5F_vs_mESC_housekeep_maternal.gr$domTSS_dist_E16_5FvsmESC[as.character(strand(double_shift_E16_5F_vs_mESC_housekeep_maternal.gr)) == "-"] <- 
    -double_shift_E16_5F_vs_mESC_housekeep_maternal.gr$domTSS_dist_E16_5FvsmESC[as.character(strand(double_shift_E16_5F_vs_mESC_housekeep_maternal.gr)) == "-"]

# single shift (maternal)
  single_shift_E16_5F_vs_mESC_housekeep.gr$domTSS_dist_E16_5FvsmESC <- single_shift_E16_5F_vs_mESC_housekeep.gr$groupX.pos - single_shift_E16_5F_vs_mESC_housekeep.gr$groupY.pos

  single_shift_E16_5F_vs_mESC_housekeep.gr$domTSS_dist_E16_5FvsmESC[as.character(strand(single_shift_E16_5F_vs_mESC_housekeep.gr)) == "-"] <- 
    -single_shift_E16_5F_vs_mESC_housekeep.gr$domTSS_dist_E16_5FvsmESC[as.character(strand(single_shift_E16_5F_vs_mESC_housekeep.gr)) == "-"]

# --- centre on mESC domTSS double shift somatic, double shift maternal an single shift granges object - Ygroup
# double shift maternal
  domTSS_mESC_double_shift_E16_5F_vs_mESC_housekeep_maternal.gr <- GRanges(seqnames = seqnames(double_shift_E16_5F_vs_mESC_housekeep_maternal.gr),
                                                                           ranges = IRanges(start = double_shift_E16_5F_vs_mESC_housekeep_maternal.gr$groupY.pos,
                                                                                            end = double_shift_E16_5F_vs_mESC_housekeep_maternal.gr$groupY.pos),
                                                                           strand = strand(double_shift_E16_5F_vs_mESC_housekeep_maternal.gr))
  mcols(domTSS_mESC_double_shift_E16_5F_vs_mESC_housekeep_maternal.gr) <- mcols(double_shift_E16_5F_vs_mESC_housekeep_maternal.gr)
  seqinfo(domTSS_mESC_double_shift_E16_5F_vs_mESC_housekeep_maternal.gr) <- seqinfo(double_shift_E16_5F_vs_mESC_housekeep_maternal.gr)
  domTSS_mESC_double_shift_E16_5F_vs_mESC_housekeep_maternal.gr$type <- "double_shift_maternal"

# double shift somatic
  domTSS_mESC_double_shift_E16_5F_vs_mESC_housekeep_somatic.gr <- GRanges(seqnames = seqnames(double_shift_E16_5F_vs_mESC_housekeep_somatic.gr),
                                                                          ranges = IRanges(start = double_shift_E16_5F_vs_mESC_housekeep_somatic.gr$groupY.pos,
                                                                                           end = double_shift_E16_5F_vs_mESC_housekeep_somatic.gr$groupY.pos),
                                                                          strand = strand(double_shift_E16_5F_vs_mESC_housekeep_somatic.gr))
  mcols(domTSS_mESC_double_shift_E16_5F_vs_mESC_housekeep_somatic.gr) <- mcols(double_shift_E16_5F_vs_mESC_housekeep_somatic.gr)
  seqinfo(domTSS_mESC_double_shift_E16_5F_vs_mESC_housekeep_somatic.gr) <- seqinfo(double_shift_E16_5F_vs_mESC_housekeep_somatic.gr)
  domTSS_mESC_double_shift_E16_5F_vs_mESC_housekeep_somatic.gr$type <- "double_shift_somatic"

# single shift (maternal)
  domTSS_mESC_single_shift_E16_5F_vs_mESC_housekeep.gr <- GRanges(seqnames = seqnames(single_shift_E16_5F_vs_mESC_housekeep.gr),
                                                                  ranges = IRanges(start = single_shift_E16_5F_vs_mESC_housekeep.gr$groupY.pos,
                                                                                   end = single_shift_E16_5F_vs_mESC_housekeep.gr$groupY.pos),
                                                                  strand = strand(single_shift_E16_5F_vs_mESC_housekeep.gr))
  mcols(domTSS_mESC_single_shift_E16_5F_vs_mESC_housekeep.gr) <- mcols(single_shift_E16_5F_vs_mESC_housekeep.gr)
  seqinfo(domTSS_mESC_single_shift_E16_5F_vs_mESC_housekeep.gr) <- seqinfo(single_shift_E16_5F_vs_mESC_housekeep.gr)
  domTSS_mESC_single_shift_E16_5F_vs_mESC_housekeep.gr$type <- "single_shift_(maternal)"

# combine all into a single GRanges
  domTSS_mESC_E16_5F_vs_mESC_all_housekeep.gr <- c(domTSS_mESC_single_shift_E16_5F_vs_mESC_housekeep.gr, domTSS_mESC_double_shift_E16_5F_vs_mESC_housekeep_somatic.gr, domTSS_mESC_double_shift_E16_5F_vs_mESC_housekeep_maternal.gr)

# convert to type to factor
  domTSS_mESC_E16_5F_vs_mESC_all_housekeep.gr$type <- factor(domTSS_mESC_E16_5F_vs_mESC_all_housekeep.gr$type, levels = unique(domTSS_mESC_E16_5F_vs_mESC_all_housekeep.gr$type))

  

# --- centre on E16.5F domTSS double shift somatic, double shift maternal an single shift granges object - Xgroup  (check code above and Y or X group)
# double shift maternal
  domTSS_E16_5F_double_shift_E16_5F_vs_mESC_housekeep_maternal.gr <- GRanges(seqnames = seqnames(double_shift_E16_5F_vs_mESC_housekeep_maternal.gr),
                                                                             ranges = IRanges(start = double_shift_E16_5F_vs_mESC_housekeep_maternal.gr$groupX.pos,
                                                                                              end = double_shift_E16_5F_vs_mESC_housekeep_maternal.gr$groupX.pos),
                                                                             strand = strand(double_shift_E16_5F_vs_mESC_housekeep_maternal.gr))
  mcols(domTSS_E16_5F_double_shift_E16_5F_vs_mESC_housekeep_maternal.gr) <- mcols(double_shift_E16_5F_vs_mESC_housekeep_maternal.gr)
  seqinfo(domTSS_E16_5F_double_shift_E16_5F_vs_mESC_housekeep_maternal.gr) <- seqinfo(double_shift_E16_5F_vs_mESC_housekeep_maternal.gr)
  domTSS_E16_5F_double_shift_E16_5F_vs_mESC_housekeep_maternal.gr$type <- "double_shift_maternal"


# double shift somatic
  domTSS_E16_5F_double_shift_E16_5F_vs_mESC_housekeep_somatic.gr <- GRanges(seqnames = seqnames(double_shift_E16_5F_vs_mESC_housekeep_somatic.gr),
                                                                            ranges = IRanges(start = double_shift_E16_5F_vs_mESC_housekeep_somatic.gr$groupX.pos,
                                                                                             end = double_shift_E16_5F_vs_mESC_housekeep_somatic.gr$groupX.pos),
                                                                            strand = strand(double_shift_E16_5F_vs_mESC_housekeep_somatic.gr))
  mcols(domTSS_E16_5F_double_shift_E16_5F_vs_mESC_housekeep_somatic.gr) <- mcols(double_shift_E16_5F_vs_mESC_housekeep_somatic.gr)
  seqinfo(domTSS_E16_5F_double_shift_E16_5F_vs_mESC_housekeep_somatic.gr) <- seqinfo(double_shift_E16_5F_vs_mESC_housekeep_somatic.gr)
  domTSS_E16_5F_double_shift_E16_5F_vs_mESC_housekeep_somatic.gr$type <- "double_shift_somatic"
  

# single shift (maternal)
  domTSS_E16_5F_single_shift_E16_5F_vs_mESC_housekeep.gr <- GRanges(seqnames = seqnames(single_shift_E16_5F_vs_mESC_housekeep.gr),
                                                                    ranges = IRanges(start = single_shift_E16_5F_vs_mESC_housekeep.gr$groupX.pos,
                                                                                     end = single_shift_E16_5F_vs_mESC_housekeep.gr$groupX.pos),
                                                                    strand = strand(single_shift_E16_5F_vs_mESC_housekeep.gr))
  mcols(domTSS_E16_5F_single_shift_E16_5F_vs_mESC_housekeep.gr) <- mcols(single_shift_E16_5F_vs_mESC_housekeep.gr)
  seqinfo(domTSS_E16_5F_single_shift_E16_5F_vs_mESC_housekeep.gr) <- seqinfo(single_shift_E16_5F_vs_mESC_housekeep.gr)
  domTSS_E16_5F_single_shift_E16_5F_vs_mESC_housekeep.gr$type <- "single_shift_(maternal)"


# combine all into a single GRanges
  domTSS_E16_5F_E16_5F_vs_mESC_all_housekeep.gr <- c(domTSS_E16_5F_single_shift_E16_5F_vs_mESC_housekeep.gr, domTSS_E16_5F_double_shift_E16_5F_vs_mESC_housekeep_somatic.gr, domTSS_E16_5F_double_shift_E16_5F_vs_mESC_housekeep_maternal.gr)

# convert to type to factor
  domTSS_E16_5F_E16_5F_vs_mESC_all_housekeep.gr$type <- factor(domTSS_E16_5F_E16_5F_vs_mESC_all_housekeep.gr$type, levels = unique(domTSS_E16_5F_E16_5F_vs_mESC_all_housekeep.gr$type))

# ---- centre to mESC domTSS and sort by distance and orientation of mESC shift 
# WW and SS heatmaps
# create windows for plotting
  up <- 500
  down <- 500
  range <- c(-up, down)
  win <- up + down

# centre on mESC domTSS
  domTSS_win.gr <- promoters(domTSS_mESC_E16_5F_vs_mESC_all_housekeep.gr, up = up, down = down)

# remove out of bound ranges
  domTSS_win.gr <- domTSS_win.gr[width(trim(domTSS_win.gr)) == win]

# randomize prior to sorting
  set.seed(10)
  random_row <- sample(length(domTSS_win.gr))
  domTSS_win.gr <- domTSS_win.gr[random_row]

# sort index - no sequences are trimmed out so I can use the .gr object for sorting - so its sorted from upstream to downstream
# first sort by distance
  sort_idx <- order(domTSS_win.gr$domTSS_dist_E16_5FvsmESC,  decreasing = FALSE)
  domTSS_win.gr <- domTSS_win.gr[sort_idx]

# second sort by type - "single_shift_(maternal)" "double_shift_somatic"    "double_shift_maternal" 
  sort_idx <- order(domTSS_win.gr$type,  decreasing = FALSE)
  domTSS_win.gr <- domTSS_win.gr[sort_idx]

# extract sequence 
  domTSS_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, domTSS_win.gr)

# plot dinucleotide plots 
  pattern <- c("WW", "SS", "TATA", "TATAWAWR")

  for (i in 1:length(pattern)) {
    hm.l = PatternHeatmap(domTSS_seq, pattern = pattern[i], 
                          coords = range, 
                          label = paste(pattern[[i]]))
    hm_smoothed.l = smoothHeatmap(hm.l, sigma = c(3, 3), output.size=c(1000, 1000))
    #scale(hm_smoothed.l) <- scale_range[[i]]
    
    pdf(paste("Figures/domTSS_mESC_GGC_DS_SS_", pattern[[i]], ".pdf", sep = "_"), height = 4, width = 4.5)
    plotHeatmapList(hm_smoothed.l, 
                    legend = FALSE, 
                    cex.label = 0.5,
                    partition = c(sum(domTSS_win.gr$type == "single_shift_(maternal)"), 
                                  sum(domTSS_win.gr$type == "double_shift_somatic"),
                                  sum(domTSS_win.gr$type == "double_shift_maternal")),
                    partition.legend = TRUE,
                    partition.lines = TRUE)
    dev.off()
  }





# import switching tag clusters
  PGC_oocyte_embryo_vs_mESC_shift.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/PGC_oocyte_embryo_vs_mESC_shift_grl.RDS")

#------------ plot CTSS coverage - E16.5 ctss coverage, centred to mESC
  CTSSnorm <- CTSSnormalizedTpm(CAGEset_PGC_embryo_merged)

# add CTSS number - 
  CTSSnorm$CTSSid <- 1:nrow(CTSSnorm)
  
  # create GRanges to overlap with consensus cluster promoters
  CTSSnorm.gr <- GRanges(seqnames = CTSSnorm$chr,
                         ranges = IRanges(start = CTSSnorm$pos,
                                          end = CTSSnorm$pos),
                         strand = CTSSnorm$strand)

# attach TPM value for E11.5 and for E16.5F
  mcols(CTSSnorm.gr) <- CTSSnorm[, "E16_5F"] 

# add seqinfo 
  seqinfo(CTSSnorm.gr) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(CTSSnorm.gr)]

# plot CTSS coverage centered on mESC, sorted by domTSS distance from mESC domTSS
# calculate weight for heatmaps coverage
  weight <-  ifelse(as.character(strand(CTSSnorm.gr)) == "+", log2(mcols(CTSSnorm.gr)[,1] + 1), -log2(mcols(CTSSnorm.gr)[,1] + 1))

# E16.5 CTSS coverage
  hm_E16_5_ctss = CoverageHeatmap(
    domTSS_win.gr,
    CTSSnorm.gr,
    coords = range,
    weight = weight,
    label = "E16_5F")

  scale(hm_E16_5_ctss) <- c(-1.5, 1.5)

# flip the negative strand to -tpm - as promoters function orients all to positive orientation, so all that maps to negative strand will be antisense 
  image(hm_E16_5_ctss)[as.character(strand(domTSS_win.gr)) == "-", ] <- -image(hm_E16_5_ctss)[as.character(strand(domTSS_win.gr)) == "-", ]

# plotting
  pdf("Figures/domTSS_coverage_mESC_GGC_DS_SS.pdf",  height = 4, width = 4.5)
  plotHeatmapList(hm_E16_5_ctss,
                  cex.label = 0.5,
                  legend = F,
                  color = c("#ca0020", "#f7f7f7", "#0571b0"),
                  partition = c(sum(domTSS_win.gr$type == "single_shift_(maternal)"), 
                                sum(domTSS_win.gr$type == "double_shift_somatic"),
                                sum(domTSS_win.gr$type == "double_shift_maternal")),
                  partition.legend = TRUE,
                  partition.lines = TRUE)
  dev.off()


# scan with nucleosome positioning pwm 
# scan the sequence with nucleosome pwm
  library(seqPattern)
  nucl_scores <- motifScanScores(domTSS_seq, motifPWM = nucleosome_pwm)  
  
  # use heatmaps for plotting
  hm.l <- new("Heatmap",
              image = nucl_scores,
              coords = as.integer(range),
              nseq = nrow(nucl_scores),
              metadata = list())
  
  hm_smoothed.l<- smoothHeatmap(hm.l, sigma = c(3, 3), output.size=c(500, 500))
  scale(hm_smoothed.l) <- c(40, 60)

# plot heatmaps
  pdf("Figures/domTSS_mESC_GGC_DS_SS_nucl_pwm.pdf", height = 4, width = 4.5)
  plotHeatmapList(hm_smoothed.l,
                  legend = T,
                  legend.pos = "l",
                  color = "Spectral",
                  cex.label = 0.5,
                  partition = c(sum(domTSS_win.gr$type == "single_shift_(maternal)"), 
                                sum(domTSS_win.gr$type == "double_shift_somatic"),
                                sum(domTSS_win.gr$type == "double_shift_maternal")),
                  partition.legend = FALSE,
                  partition.lines = TRUE)
  dev.off()


# ---- centre to E16_5 domTSS and sort by distance and orientation of mESC shift ----#
# WW and SS heatmaps
# create windows for plotting
  up <- 500
  down <- 500
  range <- c(-up, down)
  win <- up + down

# centre on mESC domTSS
  domTSS_win.gr <- promoters(domTSS_E16_5F_E16_5F_vs_mESC_all_housekeep.gr, up = up, down = down)

# remove out of bound ranges
  domTSS_win.gr <- domTSS_win.gr[width(trim(domTSS_win.gr)) == win]

# randomize prior to sorting
  set.seed(10)
  random_row <- sample(length(domTSS_win.gr))
  domTSS_win.gr <- domTSS_win.gr[random_row]

# sort index - no sequences are trimmed out so I can use the .gr object for sorting - so its sorted from upstream to downstream
# first sort by distance
  sort_idx <- order(domTSS_win.gr$domTSS_dist_E16_5FvsmESC,  decreasing = FALSE)
  domTSS_win.gr <- domTSS_win.gr[sort_idx]

# second sort by type - "single_shift_(maternal)" "double_shift_somatic"    "double_shift_maternal" 
  sort_idx <- order(domTSS_win.gr$type,  decreasing = FALSE)
  domTSS_win.gr <- domTSS_win.gr[sort_idx]

# extract sequence 
  domTSS_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, domTSS_win.gr)

# plot dinucleotide plots arround nucleosome centers (ranges are centered on nucleosome dyads)
  pattern <- c("WW", "SS", "TATA", "TATAWAWR")

  for (i in 1:length(pattern)) {
    hm.l = PatternHeatmap(domTSS_seq, pattern = pattern[i], 
                          coords = range, 
                          label = paste(pattern[[i]]))
    hm_smoothed.l = smoothHeatmap(hm.l, sigma = c(3, 3), output.size=c(1000, 1000))
    #scale(hm_smoothed.l) <- scale_range[[i]]
    
    pdf(paste("Figures/domTSS_E16_5F_GGC_DS_SS_", pattern[[i]], ".pdf", sep = "_"), height = 4, width = 4.5)
    plotHeatmapList(hm_smoothed.l, 
                    legend = FALSE, 
                    cex.label = 0.5,
                    partition = c(sum(domTSS_win.gr$type == "single_shift_(maternal)"), 
                                  sum(domTSS_win.gr$type == "double_shift_somatic"),
                                  sum(domTSS_win.gr$type == "double_shift_maternal")),
                    partition.legend = TRUE,
                    partition.lines = TRUE)
    dev.off()
  }

# import switching tag clusters
  PGC_oocyte_embryo_vs_mESC_shift.grl <- readRDS("../all_reps_PN6_2cell/intermediate_data/PGC_oocyte_embryo_vs_mESC_shift_grl.RDS")

# plot E16.5 TC coverage plot
  hm.l = CoverageHeatmap(
    domTSS_win.gr,
    track = PGC_oocyte_embryo_vs_mESC_shift.grl[["E16_5F"]],
    coords = range,
    weight = log2(PGC_oocyte_embryo_vs_mESC_shift.grl[["E16_5F"]]$groupX.tpm + 1))
  scale(hm.l) <- quantile(hm.l@image, c(0.01, 0.99))
  
  # plot heatmap
  pdf("Figures/domTSS_E16_5F_GGC_DS_SS_cov.pdf", height = 4, width = 4.5)
  plotHeatmapList(hm.l,
                  legend = FALSE,
                  color = "Greys",
                  cex.label = 0.5,
                  partition = c(sum(domTSS_win.gr$type == "single_shift_(maternal)"), 
                                sum(domTSS_win.gr$type == "double_shift_somatic"),
                                sum(domTSS_win.gr$type == "double_shift_maternal")),
                  partition.legend = TRUE,
                  partition.lines = TRUE)
  dev.off()


#------------ plot CTSS coverage - E16.5 ctss coverage, centred to mESC
CTSSnorm <- CTSSnormalizedTpm(CAGEset_PGC_embryo_merged)

# add CTSS number  
  CTSSnorm$CTSSid <- 1:nrow(CTSSnorm)

# create GRanges to overlap with consensus cluster promoters
  CTSSnorm.gr <- GRanges(seqnames = CTSSnorm$chr,
                         ranges = IRanges(start = CTSSnorm$pos,
                                          end = CTSSnorm$pos),
                         strand = CTSSnorm$strand)

# attach TPM value for E11.5 and for E16.5F
  mcols(CTSSnorm.gr) <- CTSSnorm[, "E16_5F"] 

# add seqinfo 
  seqinfo(CTSSnorm.gr) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(CTSSnorm.gr)]


# plot CTSS coverage centered on mESC, sorted by domTSS distance from mESC domTSS
# calculate weight for heatmaps coverage
  weight <-  ifelse(as.character(strand(CTSSnorm.gr)) == "+", log2(mcols(CTSSnorm.gr)[,1] + 1), -log2(mcols(CTSSnorm.gr)[,1] + 1))

# E16.5 CTSS coverage
  hm_E16_5_ctss = CoverageHeatmap(
    domTSS_win.gr,
    CTSSnorm.gr,
    coords = range,
    weight = weight,
    label = "E16_5F")
  
  scale(hm_E16_5_ctss) <- c(-1.5, 1.5)

# flip the negative strand to -tpm - as promoters function orients all to positive orientation, so all that maps to negative strand will be antisense 
  image(hm_E16_5_ctss)[as.character(strand(domTSS_win.gr)) == "-", ] <- -image(hm_E16_5_ctss)[as.character(strand(domTSS_win.gr)) == "-", ]

# plotting
  pdf("Figures/domTSS_coverage_E16_5F_GGC_DS_SS.pdf",  height = 4, width = 4.5)
  plotHeatmapList(hm_E16_5_ctss,
                  cex.label = 0.5,
                  legend = F,
                  color = c("#ca0020", "#f7f7f7", "#0571b0"),
                  partition = c(sum(domTSS_win.gr$type == "single_shift_(maternal)"), 
                                sum(domTSS_win.gr$type == "double_shift_somatic"),
                                sum(domTSS_win.gr$type == "double_shift_maternal")),
                  partition.legend = TRUE,
                  partition.lines = TRUE)
  dev.off()


# scan with nucleosome positioning pwm 
  library(seqPattern)
  nucl_scores <- motifScanScores(domTSS_seq, motifPWM = nucleosome_pwm)  

# use heatmaps for plotting
  hm.l <- new("Heatmap",
              image = nucl_scores,
              coords = as.integer(range),
              nseq = nrow(nucl_scores),
              metadata = list())

  hm_smoothed.l<- smoothHeatmap(hm.l, sigma = c(3, 3), output.size=c(500, 500))
  scale(hm_smoothed.l) <- c(40, 60)

# plot heatmaps
  pdf("Figures/domTSS_E16_5F_GGC_DS_SS_nucl_pwm.pdf", height = 4, width = 4.5)
  plotHeatmapList(hm_smoothed.l,
                  legend = T,
                  legend.pos = "l",
                  color = "Spectral",
                  cex.label = 0.5,
                  partition = c(sum(domTSS_win.gr$type == "single_shift_(maternal)"), 
                                sum(domTSS_win.gr$type == "double_shift_somatic"),
                                sum(domTSS_win.gr$type == "double_shift_maternal")),
                  partition.legend = FALSE,
                  partition.lines = TRUE)
  dev.off()
