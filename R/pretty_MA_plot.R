# Simple R script to plot an annotated MA plot 
# for a DESeq2-style DEA results table
# 
# Arguments
# results: DESeq2 results() output, can be pre-cast as data frame or not
# convert_ids: whether or not to convert from ensembl gene IDs to gene symbol
# id_col: the name of the column in your df with ensembl ids. By default, just uses rownames
# mart_name: the name of the biomart library to use. Defaults to mouse, ensembl 
# name_col: the name of the column with gene symbols to be used for plotting. Ignored if convert_ids = T
# lfc_threshold: the minimum absolute log2 fold-change for labeling. Defaults to 0 (all significant genes)
# genes: genes to subset to
# highlight_genes: gene labels to highlight with larger text
# custom_annotation: a custom gene conversion set (helpful for metagenomes); must be in biomart format

pretty_MA_plot = function(results, 
                          convert_ids = T, 
                          id_col = "row.names",
                          mart_name = "mmusculus_gene_ensembl",
                          name_col = "row.names",
                          lfc_threshold = 0,
                          genes = NULL,
                          highlight_genes = c(),
                          custom_annotation = NULL,
                          max_overlaps = 10,
                          label_alpha = 1,
                          label_text_size = (10 / .pt),
                          random_seed = 12345)
{
  require(ggplot2)
  require(ggrepel)
  require(biomaRt)
  require(dplyr)
  require(tibble)
  require(tidyverse)
  
  if(convert_ids) #from ensembl to common symbols
  {
    results = id_convert(results = results, 
               id_col = id_col,
               mart_name = mart_name,
               name_col = name_col,
               custom_annotation = custom_annotation)
    name_col = "external_gene_name"
  }
  #add colors for significant
  results = results %>% 
    mutate(color = factor(case_when(padj >= 0.05 ~ "NS",
                                    padj < 0.05 & log2FoldChange > 0 ~ "upregulated",
                                    padj < 0.05 & log2FoldChange < 0 ~ "downregulated")))
  
  plt = ggplot(results, 
               aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = color)) +
    scale_x_log10(limits = c(1, NA)) +
    scale_color_manual(values = c("NS" = alpha(colour = "grey50", alpha = 0.05),
                                  "upregulated" = "firebrick4",
                                  "downregulated" = "dodgerblue4")) +
    theme_bw(base_family = "Arial") +
    theme(legend.position = "none") + 
    xlab("Mean Expression") +
    ylab("log2(Fold Change)")
  
  if(is.null(genes))
  {
    #plot as separate up and down to prevent crossing of the origin
    plt = plt + 
      geom_label_repel(data = subset(results, padj < 0.05 & 
                                                 abs(log2FoldChange) >= lfc_threshold &
                                                 log2FoldChange > 0),
                       aes(label = .data[[name_col]], size = factor(external_gene_name %in% highlight_genes)), 
                       max.overlaps = max_overlaps,
                       fill = alpha(c("white"), label_alpha),
                       ylim = c(1, NA)) +
      geom_label_repel(data = subset(results, padj < 0.05 & 
                                       abs(log2FoldChange) >= lfc_threshold &
                                       log2FoldChange <= 0),
                       aes(label = .data[[name_col]], size = factor(external_gene_name %in% highlight_genes)), 
                       max.overlaps = max_overlaps, 
                       fill = alpha(c("white"), label_alpha),
                       ylim = c(NA, -1)) +
      scale_size_manual(values = c("TRUE" = label_text_size * 1.5,
                                   "FALSE" = label_text_size))
  } else
  {
    plt = plt + 
      geom_label_repel(data = subset(results, external_gene_name %in% genes & log2FoldChange > 0),
                       aes(label = .data[[name_col]], size = factor(external_gene_name %in% highlight_genes)),
                       min.segment.length = 0.1, 
                       max.overlaps = max_overlaps, 
                       fill = alpha(c("white"), label_alpha),
                       ylim = c(1, NA)) +
      geom_label_repel(data = subset(results, external_gene_name %in% genes & log2FoldChange <= 0),
                       aes(label = .data[[name_col]], size = factor(external_gene_name %in% highlight_genes)),
                       min.segment.length = 0.1, 
                       max.overlaps = max_overlaps, 
                       fill = alpha(c("white"), label_alpha),
                       ylim = c(NA, -1)) +
      scale_size_manual(values = c("TRUE" = label_text_size * 1.5,
                                     "FALSE" = label_text_size))
  }
                                 
  
  return(plt)
}   

id_convert = function(results,
                      id_col = "row.names",
                      mart_name = "mmusculus_gene_ensembl",
                      name_col = "row.names",
                      custom_annotation = NULL)
{
  require(biomaRt)
  if(!is.data.frame(results)) #convert to standard data frame, as necessary
  {
    results = as.data.frame(results)
  }
  
  if(!is.null(custom_annotation))
  {
    conv = custom_annotation
  } else
  {
    mart = useMart("ensembl", mart_name)
    conv = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                 mart = mart)
  }
  results = rownames_to_column(results, var = id_col)
  results = left_join(results,
                      conv,
                      by = setNames("ensembl_gene_id", id_col))
  return(results)
}