#' Function to make all pairwise comparisons for DESeq2 DEA
#' 
#' @param des DESeq2 object with or without DESeq() function already performed (will be rerun either way).
#' @param change_design if true, design(des) is changed to ~comparison_col. FALSE will use the current design.
#' @param comparison_col column to use for DEA. Defaults to as.character(design(des))[2].
#' @param fit_type DESeq2 dispersion fitting model
#' @param min_reps_replacement minimum number of replicates for gene replacement in DESeq2
#' @param sort_direction sort in ascending or descending alphanumeric order? 
#' @param df_output if true, convert DESeq results objects to data frames
#' @param save_csv if true, output individual CSV for each comparison to output_directory
#' @param save_pdf if true, output PDFs of MA plots to output_directory
#' @param output_directory full path to where to save files if output is enabled. Created if it does not already exist. If NA, no output.
#' @param cores number of cores to use for DEA
#' @param convert_ids for pretty_ma_plot
#' @param id_col for pretty_ma_plot
#' @param mart_name for pretty_ma_plot
#' @param name_col for pretty_ma_plot
#' @param lfc_threshold for pretty_ma_plot
#' @param genes for pretty_ma_plot
#' @param custom_annotation for pretty_ma_plot
#' @param max_overlaps for pretty_ma_plot
#' @param label_alpha for pretty_ma_plot
#' @param label_text_size for pretty_ma_plot
#' @param random_seed for pretty_ma_plot
#' @param pdf_width width of pdf in inches
#' @param pdf_height height of pdf in inches
#' @param y_min minimum value of y for resultant plots
#' @param y_max maximum value of y for resultant plots
#' @param label_only_sig if true, label only significant genes from 'genes' argument (passed to pretty_MA_plot)
#' @param highlight_genes gene labels to highlight with larger text (passed to pretty_MA_plot)
#' @param label_oor if true, adds triangles representing genes out of the plot limits. Defaults to FALSE. (passed to pretty_MA_plot)
#' @param alpha value of alpha for differential expression analysis. Defaults to 0.05
#' @return a list of lists. "MA" is a list of ggplot2-editable MA plots. "hits" is a list of results objects.
#' @import DESeq2 BiocParallel dplyr magrittr biomaRt tibble
#' @export

get_pairwise_DESeq = function(des, change_design = TRUE, comparison_col = NA, fit_type = "parametric", 
                              min_reps_replacement = 7, sort_direction = c("factor_order", "ascending", "descending"),
                              save_csv = F, save_pdf = F, df_output = F, 
                              output_directory = NA, cores = 1, convert_ids = T, 
                              id_col = "row.names", mart_name = "mmusculus_gene_ensembl",
                              name_col = "row.names", lfc_threshold = 0, genes = NULL,
                              custom_annotation = NULL, max_overlaps = 10, label_alpha = 1,
                              random_seed = 12345, pdf_width = 6, pdf_height = 4, label_text_size = (10 / .pt),
                              y_min = NA, y_max = NA, label_only_sig = FALSE, label_oor = FALSE,
                              highlight_genes = c(), alpha = 0.05)
{
  
  if(!("DESeq2" %in% .packages()))
  {
    library(DESeq2)
  }
  if(!("BiocParallel" %in% .packages()))
  {
    library(BiocParallel)
  }
  if(!("dplyr" %in% .packages()))
  {
    library(dplyr)
  }
  if(!("magrittr" %in% .packages()))
  {
    library(magrittr)
  }
  if(!("biomaRt" %in% .packages()))
  {
    library(biomaRt)
  }
  if(!("tibble" %in% .packages()))
  {
    library(tibble)
  }
  
  register(MulticoreParam(cores))
  
  #evaluate sort direction
  sort_direction = match.arg(sort_direction)
  
  #assume one-factor design if not specified
  if(is.na(comparison_col))
  {
    comparison_col = as.character(design(des))[2]
  }
  
  #if user wants design change
  if(change_design == TRUE)
  {
    design(des) = as.formula(paste0("~", comparison_col))
  }
  
  #run DEA
  dge = DESeq(des, fitType = fit_type, minReplicatesForReplace = min_reps_replacement, parallel = T)
  
  #get all possible combinations and sort
  all_comparisons = combn(levels(dge[[comparison_col]]), 2, simplify = F)
  all_comparisons = lapply(all_comparisons, function(comp){
    if(sort_direction == "factor_order")
    {
      comp = factor(comp, levels = levels(dge[[comparison_col]]))
      comp = sort(comp)
    }
    else if(sort_direction == "ascending")
    {
      comp = sort(comp)
    } else if(sort_direction == "descending")
    {
      comp = sort(comp, decreasing = T)
    } else
    {
      stop("Error: invalid sort order selection")
    }
    comp = c(comparison_col, as.character(comp)) #get in format for results()
    return(comp) })
  
  #perform comparisons
  res = lapply(all_comparisons, function(comp){
    hits = results(dge, contrast = comp, alpha = alpha, parallel = T)
    return(hits) })
  names(res) = vapply(all_comparisons, function(comp){
    name = paste(comp[1], comp[2], "over", comp[3], sep = "_")
    return(name) }, FUN.VALUE = "char")
  
  #make MA plots
  ma_plots = lapply(res, function(comp){
    plot = pretty_MA_plot(results = comp, convert_ids = convert_ids, 
                          id_col = id_col, mart_name = mart_name,
                          name_col = name_col, lfc_threshold = lfc_threshold, genes = genes,
                          custom_annotation = custom_annotation, max_overlaps = max_overlaps, 
                          label_alpha = label_alpha, random_seed = random_seed, label_text_size = label_text_size,
                          y_min = y_min, y_max = y_max, label_only_sig = label_only_sig,
                          label_oor = label_oor, highlight_genes = highlight_genes)
    return(plot) })
  names(ma_plots) = names(res)
  
  #convert results to dataframes for easier handling
  if(!is.null(custom_annotation))
  {
    gene_conv = custom_annotation
  } else
  {
    mart = useMart("ensembl", mart_name)
    gene_conv = getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = mart)
  }
    
  res_df = lapply(res, function(comp){
    df = as.data.frame(comp) %>% 
      rownames_to_column("ensembl_gene_id") %>% 
      left_join(., gene_conv) %>% 
      dplyr::relocate(ensembl_gene_id, external_gene_name)
    return(df) })
  
  #data output   
  if(!is.na(output_directory))
  {
    if(!dir.exists(output_directory))
    {
      dir.create(output_directory)
    }
    
    if(save_csv == T)
    {
      for(i in 1:length(res_df))
      {
        outpath = paste0(output_directory, "/", names(res_df)[i], ".csv")
        write.csv(res_df[[i]], outpath)
      }
    }
    if(save_pdf == T)
    {
      for(i in 1:length(ma_plots))
      {
        outpath = paste0(output_directory, "/", names(ma_plots)[i], ".pdf")
        cairo_pdf(outpath,
                 width = pdf_width,
                 height = pdf_height,
                 family = "Arial")
        plot(ma_plots[[i]])
        dev.off()
      }
    }
  }
    
  #return data
  if(df_output == T)
  {
    output = list("hits" = res_df, "MA" = ma_plots)
  } else
  {
    output = list("hits" = res, "MA" = ma_plots)
  }
  return(output)
}