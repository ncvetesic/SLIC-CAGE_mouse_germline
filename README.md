Global regulatory transitions at core promoters demarcate the mammalian
germline cycle
================
Nevena Cvetesic
06/11/2020

## Abstract

Core promoters integrate regulatory inputs of genes. Global dynamics of promoter usage can reveal systemic changes in how genomic sequence is interpreted by the cell. Here we report the first analysis of promoter dynamics and code switching in the mammalian germ line, characterising the full cycle of transitions from embryonic stem cells through germline, oogenesis, and zygotic genome activation. Using Super Low Input Carrier-CAGE (SLIC-CAGE) we show that mouse germline development starts with the somatic promoter code, followed by a prominent switch to the maternal code during follicular oogenesis. The sequence features underlying the shift from somatic to maternal code are conserved across vertebrates, despite large differences in promoter nucleotide compositions. In addition, we show that, prior to this major shift, the promoters of gonadal germ cells diverge from the canonical somatic transcription initiation. This divergence is distinct from the promoter code used later by developing oocytes and reveals genome-wide promoter remodelling associated with alternative nucleosome positioning during early female and male germline development. Collectively, our findings establish promoter-level regulatory transitions as a central, conserved feature of the vertebrate life cycle.

## Package requirements

  - Install the following R packages from CRAN:   
  dplyr\_0.8.4, RColorBrewer\_1.1-2, ggseqlogo\_0.1, ggplot2\_3.2.1, corrplot\_0.84

  - Install the following R packages from Bioconductor:  
  CAGEr\_1.20.0, BSgenome\_1.52.0, BSgenome.Mmusculus.UCSC.mm10\_1.4.0,
    rtracklayer\_1.44.4, Biostrings\_2.52.0, org.Mm.eg.db\_3.8.2,
    TxDb.Mmusculus.UCSC.mm10.knownGene\_3.4.7,GenomicFeatures\_1.36.4,
    AnnotationDbi\_1.46.1, Biobase\_2.44.0, ChIPseeker\_1.20.0,
    TFBSTools\_1.22.0, heatmaps\_1.8.0, Gviz\_1.29.1,
    GenomicRanges\_1.36.1, GenomeInfoDb\_1.20.0, IRanges\_2.18.3,
    BiocGenerics\_0.30.0, DESeq2\_1.24.0, seqPattern\_1.16.0

## Figure to code map

  - Figure 1[C](analysis/01_CTSS_expression_correlation.R)
  - Figure 1[D](analysis/02_CTSS_PCA.R)
  - Figure 2[B](analysis/03_domTSS_distr_distribution.R)
  - Figure 2[C-F](analysis/04_heatmaps.R), 3[B,C](analysis/04_heatmaps.R), 4[D,F,G](analysis/04_heatmaps.R)
  - Figure 2[G](analysis/05_TBPpwm_match_distribution.R), 3[A,E](analysis/05_TBPpwm_match_distribution.R)  
  - Figure 2[H](analysis/06_seqlogos.R), 3[D](analysis/06_seqlogos.R)
  - Figure 2[J](analysis/07_Wbox_stretch_length.R)
  - Figure 4[E](analysis/08_IQwidth_correlation.R)
  
## Other code  

  - identifying [shifting promoters](analysis/shifting_promoters.R)
  - zebrafish [shifting promoters](analysis/zebrafish_CAGE_shifts.R)

