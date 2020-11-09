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
# custom_annotation: a custom gene conversion set (helpful for metagenomes); must be in biomart format

pretty_MA_plot = function(results, 
                          convert_ids = T, 
                          id_col = "row.names",
                          mart_name = "mmusculus_gene_ensembl",
                          name_col = "row.names",
                          lfc_threshold = 0,
                          genes = NULL,
                          custom_annotation = NULL)
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
                                    padj < 0.05 & log2FoldChange >= lfc_threshold ~ "upregulated",
                                    padj < 0.05 & log2FoldChange <= -lfc_threshold ~ "downregulated")))
  
  plt = ggplot(results, 
               aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = color)) +
    scale_x_log10(limits = c(1, NA)) +
    scale_color_manual(values = c("NS" = alpha(colour = "grey50", alpha = 0.05),
                                  "upregulated" = "firebrick4",
                                  "downregulated" = "dodgerblue4")) +
    theme(legend.position = "none") + 
    xlab("Mean Expression") +
    ylab("log2(Fold Change)")
  
  if(is.null(genes))
  {
    plt = plt + geom_label_repel(data = subset(results, padj < 0.05 & abs(log2FoldChange) >= lfc_threshold),
                     aes_string(label = name_col))
  } else
  {
    plt = plt + geom_label_repel(data = subset(results, external_gene_name %in% genes),
                                 aes_string(label = name_col),
                                 min.segment.length = 0.1)
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
  
  mart = useMart("ensembl", mart_name)
  if(!is.null(custom_annotation))
  {
    conv = custom_annotation
  } else
  {
    conv = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                 mart = mart)
  }
  results = rownames_to_column(results, var = id_col)
  results = left_join(results,
                      conv,
                      by = setNames("ensembl_gene_id", id_col))
  return(results)
}