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

pretty_MA_plot = function(results, 
                          convert_ids = T, 
                          id_col = "row.names",
                          mart_name = "mmusculus_gene_ensembl",
                          name_col = "row.names",
                          lfc_threshold = 0)
{
  require(ggplot2)
  require(ggrepel)
  require(biomaRt)
  require(dplyr)
  require(tibble)
  
  if(convert_ids) #from ensembl to common symbols
  {
    results = id_convert(results = results, 
               id_col = id_col,
               mart_name = mart_name,
               name_col = name_col)
    name_col = "external_gene_name"
  }
  
  plt = ggplot(results, 
               aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) >= lfc_threshold)) +
    geom_label_repel(data = subset(results, padj < 0.05 & abs(log2FoldChange) >= lfc_threshold),
                     aes_string(label = name_col)) +
    scale_x_log10(limits = c(1, NA)) +
    scale_color_manual(values = c("grey50", "firebrick4")) +
    theme(legend.position = "none") + 
    xlab("Mean Expression") +
    ylab("log2(Fold Change)")
  
  suppressWarnings(print(plt))
}   

id_convert = function(results,
                      id_col = "row.names",
                      mart_name = "mmusculus_gene_ensembl",
                      name_col = "row.names")
{
  if(!is.data.frame(results)) #convert to standard data frame, as necessary
  {
    results = as.data.frame(results)
  }
  
  mart = useMart("ensembl", mart_name)
  conv = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
               mart = mart)
  results = rownames_to_column(results, var = id_col)
  results = left_join(results,
                      conv,
                      by = setNames("ensembl_gene_id", id_col))
  return(results)
}