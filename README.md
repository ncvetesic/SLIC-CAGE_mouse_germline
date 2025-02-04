Global regulatory transitions at core promoters demarcate the mammalian
germline cycle
================
Nevena Cvetesic
06/11/2020

## Abstract

Core promoters integrate regulatory inputs of genes. Global dynamics of promoter usage can reveal systemic changes in how genomic sequence is interpreted by the cell. Here we report the first analysis of promoter dynamics and code switching in the mammalian germ line, characterising the full cycle of transitions from embryonic stem cells through germline, oogenesis, and zygotic genome activation. Using Super Low Input Carrier-CAGE (SLIC-CAGE) we show that mouse germline development starts with the somatic promoter code, followed by a prominent switch to the maternal code during follicular oogenesis. The sequence features underlying the shift from somatic to maternal code are conserved across vertebrates, despite large differences in promoter nucleotide compositions. In addition, we show that, prior to this major shift, the promoters of gonadal germ cells diverge from the canonical somatic transcription initiation. This divergence is distinct from the promoter code used later by developing oocytes and reveals genome-wide promoter remodelling associated with alternative nucleosome positioning during early female and male germline development. Collectively, our findings establish promoter-level regulatory transitions as a central, conserved feature of the vertebrate life cycle.

## Package requirements

  - Install the following R packages from CRAN:   
    - dplyr_0.8.4
    - RColorBrewer_1.1-2
    - ggseqlogo_0.1 
    - ggplot2_3.2.1
    - corrplot_0.84
    - seqinr_3.6-1
    - matrixStats_0.55.0
    - kohonen_3.0.10
    - stringr_1.4.0
    - forcats_0.5.0 
    - cowplot_1.0.0
    - magrittr_1.5 
    - readxl_1.3.1 
    
  - Install the following R packages from Bioconductor:  
    - CAGEr_1.20.0
    - BSgenome_1.52.0
    - BSgenome.Mmusculus.UCSC.mm10_1.4.0,
    - rtracklayer_1.44.4
    - Biostrings_2.52.0
    - org.Mm.eg.db_3.8.2,
    - TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.7
    - GenomicFeatures_1.36.4,
    - AnnotationDbi_1.46.1
    - Biobase_2.44.0 
    - ChIPseeker_1.20.0,
    - TFBSTools_1.22.0
    - heatmaps_1.8.0
    - Gviz_1.29.1,
    - GenomicRanges_1.36.1
    - GenomeInfoDb_1.20.0
    - IRanges_2.18.3,
    - BiocGenerics_0.30.0 
    - DESeq2_1.24.0
    - seqPattern_1.16.0
    - ChIPpeakAnno_3.18.2
    - clusterProfiler_3.12.0
    - Biostrings_2.52.0
    
    
## Figure to code map
### Main Figures

  - Figure 1[C](analysis/01_CTSS_expression_correlation.R)
  - Figure 1[D](analysis/02_TC_tpm_PCA.R)
  - Figure 2[B](analysis/03_domTSS_dist_distribution.R)
  - Figure 2[C-F](analysis/04_heatmaps.R), 3[B,C](analysis/04_heatmaps.R), 4[D,F,G](analysis/04_heatmaps.R)
  - Figure 2[G](analysis/05_TBPpwm_match_distribution.R), 3[A,E](analysis/05_TBPpwm_match_distribution.R)  
  - Figure 2[H](analysis/06_seqlogos.R), 3[D](analysis/06_seqlogos.R)
  - Figure 2[J](analysis/07_Wbox_stretch_length.R)
  - Figure 4[E](analysis/08_IQwidth_correlation.R)

### Extended Data Figures
  - Extended Data Figure 2[A](analysis/09_IQwidth_distribution.R)
  - Extended Data Figure 2[B](analysis/10_IQwidth_distribution_boxplot.R)
  - Extended Data Figure 2[C](analysis/11_narrow_broad_promoters_no.R)
  - Extended Data Figure 2[D](analysis/12_promoter_genomic_locations.R)
  - Extended Data Figure 2[E,F](analysis/13_CTSS_PCA.R)
  - Extended Data Figure 3[A-D](analysis/14_distance_annotation.R)
  - Extended Data Figure 4[A,B](analysis/15_GO_enrichment.R), 12[A,B](analysis/15_GO_enrichment.R)
  - Extended Data Figure 4[C](analysis/16_GC_promoter_content.R), 12[C](analysis/16_GC_promoter_content.R)
  - Extended Data Figure 4[D](analysis/17_tpm_distribution_shifts.R), 12[D](analysis/17_tpm_distribution_shifts.R)
  - Extended Data Figure 4[E](analysis/18_Venn_shifts_overlap.R), 12[E,F](analysis/18_Venn_shifts_overlap.R)
  - Extended Data Figure 4[F](analysis/03_domTSS_dist_distribution.R)
  - Extended Data Figure 5[B,D,E](analysis/04_heatmaps.R), 6[A-D](analysis/04_heatmaps.R)
  - Extended Data Figure 7[A-D](analysis/19_tetranucleotide_composition.R)
  - Extended Data Figure 8[A](analysis/20_SOM_promoter_classification.R)
  - Extended Data Figure 8[B](analysis/21_SOM_promoter_class_genomic_location.R)
  - Extended Data Figure 8[C](analysis/22_SOM_class_expression_correlation.R)
  - Extended Data Figure 8[D-G](analysis/23_SOM_class_GO_enrichment.R)
  - Extended Data Figure 9[A,C](analysis/05_TBPpwm_match_distribution.R)
  - Extended Data Figure 9[B](analysis/22_SOM_class_expression_correlation.R)
  - Extended Data Figure 9[D,F,G](analysis/04_heatmaps.R)
  - Extended Data Figure 9[E](analysis/06_seqlogos.R)
  - Extended Data Figure 10[A,B](analysis/04_heatmaps.R)
  - Extended Data Figure 11[A,B](analysis/24_somaticTSS_expression.R)
  - Extended Data Figure 12[G](analysis/25_GGC_shifts_expression_corr.R)
  - Extended Data Figure 12[H](analysis/26_CpG_correlation_shift_spread.R)
  - Extended Data Figure 12[I](analysis/27_GGC_spread_methylation.R)
  - Extended Data Figure 12[J,K](analysis/28_WW_periodicity_GGC_shifts.R)
  - Extended Data Figure 12[L](analysis/29_WW_autocorrelation_GGC_shifts.R)
  - Extended Data Figure 12[M](analysis/30_WW_autocorrelation_mESC.R)
  - Extended Data Figure 12[N](analysis/04_heatmaps.R)
  - Extended Data Figure 13[A-D](analysis/31_pre_nucleosome_pwm.R)
  - Extended Data Figure 13[E](analysis/32_GGC_shift_correlations.R)
  - Extended Data Figure 14[A](analysis/33_GGC_single_double_shift_stratification.R)
  - Extended Data Figure 14[E,F](analysis/34_GGC_single_double_shift_heatmaps.R)
  - Extended Data Figure 15[A](analysis/35_SOM_domCTSS_classification.R)
  - Extended Data Figure 15[B,C](analysis/36_SOM_domCTSS_WW_periodicity.R)
  - Extended Data Figure 15[D,E](analysis/37_SOM_domCTSS_WW_autocorrelation.R)
  - Extended Data Figure 16[A](analysis/38_Tbpl2_expression.R)
  - Extended Data Figure 16[B](analysis/39_H3f3ab_expression.R)
  - Extended Data Figure 16[C](analysis/40_H3_3_histone_KO.R)
  - Extended Data Figure 16[D](analysis/41_altGTF_expression.R)
  
## Other code  
  - general processing of CAGE libraries starting from [bam files](analysis/CAGE_processing.R)
  - identifying [shifting promoters](analysis/shifting_promoters.R)
  - zebrafish [shifting promoters](analysis/zebrafish_CAGE_shifts.R)

