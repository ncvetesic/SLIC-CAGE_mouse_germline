#/*==========================================================================#*/
#' ## Extended Data Figure 16A
#/*==========================================================================#*/
# Tbpl2 expression

# load ta cluster data
  tc_PGCs_oocyte_embryo_anno.grl <- readRDS( "../all_reps_PN6_2cell/intermediate_data/tc_PGCs_oocyte_embryo_anno_grl.RDS")

# select for the transcript
  tbpl2.grl <- sapply(tc_PGCs_oocyte_embryo_anno.grl, function(x) x[x$transcriptId == "ENSMUST00000080453.7"])

# exclude TBP2KO
  tbpl2.grl <- tbpl2.grl[-c(15, 17)]
  tbpl2.df <- data.frame( samples = names(tbpl2.grl) , tpm =  sapply(tbpl2.grl, function(x) return(sum(x$tpm))))
  tbpl2.df$samples <- factor(tbpl2.df$samples, levels = rev(names(tbpl2.grl)))

# plot
  library(ggplot2)
  pA <- ggplot(data = tbpl2.df, aes(x = samples,
                                    y = tpm,
                                    fill = samples)) + 
    geom_bar(stat = "identity", colour = "black", alpha = 0.6) + 
    scale_fill_manual(values = rev(c("firebrick1",  "#FDE725FF",  "#E8E419FF",  "#D1E11CFF",
                                     "#B9DE28FF", "#A2DA37FF",  "#8BD646FF",  "#75D054FF",
                                     "#61CB5FFF", "#4FC46AFF" ,"#3FBC73FF","#31B57BFF",
                                     "#440154FF", "#481B6DFF", "#46337EFF","#3F4889FF",
                                     "#365C8DFF", "#2E6E8EFF"))) +
    theme_light() +
    theme(text = element_text(size = 20, colour = "black"), 
          axis.text.x = element_text(size = 20, colour = "black"),
          axis.text.y = element_text(size = 20, colour = "black"),
          legend.position = "none") +
    coord_flip()
  pdf("Figures/TBP2_expression.pdf", width = 6, height = 6)
  print(pA)
  dev.off() 