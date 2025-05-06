#' Simple R script to plot an annotated MA plot for a DESeq2-style DEA results table
#' 
#' @param results DESeq2 results() output, can be pre-cast as data frame or not
#' @param convert_ids whether or not to convert from ensembl gene IDs to gene symbol
#' @param id_col the name of the column in your df with ensembl ids. By default, just uses rownames
#' @param mart_name the name of the biomart library to use. Defaults to mouse, ensembl 
#' @param name_col the name of the column with gene symbols to be used for plotting. Ignored if convert_ids = T
#' @param lfc_threshold the minimum absolute log2 fold-change for labeling. Defaults to 0 (all significant genes)
#' @param genes genes to subset to
#' @param highlight_genes gene labels to highlight with larger text
#' @param custom_annotation a custom gene conversion set (helpful for metagenomes); must be in biomart format
#' @param max_overlaps passed to geom_label_repel()
#' @param random_seed passed to geom_label_repel
#' @param label_alpha alpha value for gene labels
#' @param label_text_size text size of gene labels
#' @param y_min minimum value of y for resultant plot
#' @param y_max maximum value of y for resultant plot
#' @param label_only_sig if true, label only significant genes from 'genes' argument
#' @param label_oor if true, adds triangles representing genes out of the plot limits. Defaults to FALSE.
#' @param alpha value of alpha for differential expression analysis. Defaults to 0.05
#' @import DESeq2 ggplot2 ggrepel biomaRt
#' @return an MA plot generated in ggplot2 dplyr magrittr tibble
#' @export

pretty_MA_plot = function(results, 
                          convert_ids = TRUE, 
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
                          y_min = NA,
                          y_max = NA,
                          label_only_sig = FALSE,
                          label_oor = FALSE,
                          random_seed = 12345,
                          alpha = 0.05)
{
  
  if(!("DESeq2" %in% .packages()))
  {
    library(DESeq2)
  }
  if(!("dplyr" %in% .packages()))
  {
    library(dplyr)
  }
  if(!("magrittr" %in% .packages()))
  {
    library(magrittr)
  }
  if(!("tibble" %in% .packages()))
  {
    library(tibble)
  }
  if(!("ggplot2" %in% .packages()))
  {
    library(ggplot2)
  }
  if(!("ggrepel" %in% .packages()))
  {
    library(ggrepel)
  }
  if(!("biomaRt" %in% .packages()))
  {
    library(biomaRt)
  }
  
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
    mutate(color = factor(case_when(padj >= alpha ~ "NS",
                                    padj < alpha & log2FoldChange > 0 ~ "upregulated",
                                    padj < alpha & log2FoldChange < 0 ~ "downregulated")),
           #effectively ignored if label_oor == FALSE, but label for later shape change
           oor = factor(case_when(log2FoldChange > y_max ~ "OOR High",
                                  log2FoldChange < y_min ~ "OOR Low",
                                  TRUE ~ "In Range")))
  
  #set to outliers to min/max yval
  if(label_oor == TRUE)
  {
    results = results %>% 
      dplyr::mutate(log2FoldChange = case_when(oor == "OOR High" ~ y_max,
                                               oor == "OOR Low" ~ y_min,
                                               oor == "In Range" ~ log2FoldChange))
  }
  
  #split highlight genes into up-, downregulated
  highlight_genes_up = 
  
  plt = ggplot(results, 
               aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = color, shape = oor)) +
    scale_x_log10(limits = c(1, NA)) +
    scale_color_manual(values = c("NS" = alpha(colour = "grey50", alpha = 0.05),
                                  "upregulated" = "firebrick4",
                                  "downregulated" = "dodgerblue4")) +
    scale_shape_manual(values = c("In Range" = 19,
                                  "OOR High" = 2,
                                  "OOR Low" = 6)) +
    theme_bw(base_family = "Arial") +
    theme(legend.position = "none") + 
    xlab("Mean Expression") +
    ylab("log2(Fold Change)")
  
  if(is.null(genes))
  {
    #plot as separate up and down to prevent crossing of the origin
    plt = plt + 
      #upregulated
      geom_label_repel(data = subset(results, padj < alpha & 
                                       abs(log2FoldChange) >= lfc_threshold &
                                       log2FoldChange > 0),
                       aes(label = .data[[name_col]]), 
                       size = label_text_size,
                       max.overlaps = max_overlaps,
                       fill = alpha(c("white"), label_alpha),
                       ylim = c(1, NA)) +
      #downregulated
      geom_label_repel(data = subset(results, padj < alpha & 
                                       abs(log2FoldChange) >= lfc_threshold &
                                       log2FoldChange <= 0),
                       aes(label = .data[[name_col]]), 
                       size = label_text_size,
                       max.overlaps = max_overlaps, 
                       fill = alpha(c("white"), label_alpha),
                       ylim = c(NA, -1)) +
      #highlight genes (priority): up
      geom_label_repel(data = subset(results, external_gene_name %in% highlight_genes & log2FoldChange > 0),
                       aes(label = .data[[name_col]]),
                       size = label_text_size * 1.5,
                       min.segment.length = 0.1, 
                       max.overlaps = 1e3, #arbitrarily high
                       fill = alpha(c("white"), label_alpha),
                       ylim = c(1, NA)) +
    #highlight genes (priority): down
    geom_label_repel(data = subset(results, external_gene_name %in% highlight_genes & log2FoldChange < 0),
                     aes(label = .data[[name_col]]),
                     size = label_text_size * 1.5,
                     min.segment.length = 0.1, 
                     max.overlaps = 1e3, #arbitrarily high
                     fill = alpha(c("white"), label_alpha),
                     ylim = c(NA, -1)) 
  } else
  {
    if(label_only_sig == TRUE)
    {
      hits = results %>% 
        dplyr::filter(padj < alpha) %>% 
        .$external_gene_name
      genes = intersect(genes, hits)
      highlight_genes = intersect(highlight_genes, hits)
    }
    plt = plt + 
      #upregulated
      geom_label_repel(data = subset(results, external_gene_name %in% genes & log2FoldChange > 0),
                       aes(label = .data[[name_col]]),
                       size = label_text_size,
                       min.segment.length = 0.1, 
                       max.overlaps = max_overlaps, 
                       fill = alpha(c("white"), label_alpha),
                       ylim = c(1, NA)) +
      #downregulated
      geom_label_repel(data = subset(results, external_gene_name %in% genes & log2FoldChange <= 0),
                       aes(label = .data[[name_col]]),
                       size = label_text_size,
                       min.segment.length = 0.1, 
                       max.overlaps = max_overlaps, 
                       fill = alpha(c("white"), label_alpha),
                       ylim = c(NA, -1)) +
      #highlight genes (priority): up
      geom_label_repel(data = subset(results, external_gene_name %in% highlight_genes & log2FoldChange > 0),
                       aes(label = .data[[name_col]]),
                       size = label_text_size * 1.5,
                       min.segment.length = 0.1, 
                       max.overlaps = 1e3, #arbitrarily high
                       fill = alpha(c("white"), label_alpha),
                       ylim = c(1, NA)) +
      geom_label_repel(data = subset(results, external_gene_name %in% highlight_genes & log2FoldChange < 0),
                       aes(label = .data[[name_col]]),
                       size = label_text_size * 1.5,
                       min.segment.length = 0.1, 
                       max.overlaps = 1e3, #arbitrarily high
                       fill = alpha(c("white"), label_alpha),
                       ylim = c(NA, -1)) 
  }
  
  #edit ylim as necessary
  if(!is.na(y_min) || !is.na(y_max))
  {
    plt = plt +
      ylim(y_min, y_max)
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