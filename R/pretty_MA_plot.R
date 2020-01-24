# Simple R script to plot an annotated MA plot 
# for a DESeq2-style DEA results table
# 
# Arguments
# results_df: DESeq2 results() output, can be pre-cast as data frame or not
# convert_ids: whether or not to convert from ensembl gene IDs to gene symbol
# id_col: the name of the column in your df with ensembl ids. By default, just uses rownames
# mart_name: the name of the biomart library to use. Defaults to mouse, ensembl 
# name_col: the name of the column with gene symbols to be used for plotting. Ignored if convert_ids = T

pretty_MA_plot = function(results_df, 
                          convert_ids = T, 
                          id_col = "row.names",
                          mart_name = "mmusculus_gene_ensembl",
                          name_col = "row.names")
{
  require(ggplot2)
  require(ggrepel)
  require(biomaRt)
  
  if(!is.data.frame(results_df)) #convert to standard data frame, as necessary
  {
    results_df = as.data.frame(results_df)
  }
  
  if(convert_ids) #from ensembl to common symbols
  {
    mart = useMart("ensembl", mart_name)
    conv = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                 mart = mart)
    results_df = merge(results_df,
                       conv,
                       all.x = T,
                       all.y = F,
                       by.x = id_col,
                       by.y = "ensembl_gene_id")
    name_col = "ensembl_gene_id"
  }
  
  plt = ggplot(results_df, 
               aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = padj < 0.05)) +
    geom_label_repel(data = subset(results_df, padj < 0.05),
                     aes_string(label = name_col)) +
    scale_x_log10(limits = c(1, NA)) +
    scale_color_manual(values = c("grey50", "firebrick4")) +
    theme(legend.position = "none") + 
    xlab("log10(Mean Expression)") +
    ylab("log2(Fold Change)")
  
  print(plt)
}