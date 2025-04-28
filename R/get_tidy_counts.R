#' script to extract counts for specified genes from a DESeq object, with metadata
#' 
#' @param obj DESeq2 object
#' @param goi vector of genes of interest in ensembl or gene name format; if NULL (default) returns all genes
#' @param goi_format "ensembl_gene_id" (default) or "external_gene_name" (symbols)
#' @param species species in ensembl format e.g. mmusculus (default); ignored if custom annotation is supplied
#' @param custom_annotation custom conversion between ensembl and external gene name
#' @param df_shape shape of output data frame. "long": unique row for each gene * sample combination. "wide": each gene (symbol) gets a column
#' @import DESeq2 dplyr magrittr tibble
#' @return a dataframe of counts for specified genes from a DESeq object, with metadata included
#' @export

get_tidy_counts = function(obj,
                           goi = NULL,
                           goi_format = "ensembl_gene_id",
                           species = "mmusculus",
                           custom_annotation = NULL,
                           df_shape = "long")
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
    library(maggritr)
  }
  if(!("tibble" %in% .packages()))
  {
    library(tibble)
  }
  
  if(is.null(custom_annotation))
  {
    require(biomaRt)
    mart_name = paste0(species, "_gene_ensembl")
    mart = useMart("ensembl", mart_name)
    gene_conv = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      mart = mart) %>% 
      dplyr::filter(ensembl_gene_id %in% rownames(obj))
  } else
  {
    gene_conv = custom_annotation
  }
  
  if(!is.null(goi))
  {
    if(goi_format == "external_gene_name")
    {
      ensembl_goi = gene_conv %>% 
        dplyr::filter(external_gene_name %in% goi) %>% 
        .$ensembl_gene_id
    } else if(goi_format == "ensembl_gene_id")
    {
      ensembl_goi = goi
    } else
    {
      stop("Error: goi must be in ensembl (ensembl_gene_id) or gene symbol (external_gene_id) format")
    }
  } else #in this case we take all genes
  {
    ensembl_goi = gene_conv$ensembl_gene_id
  }
  
  genes = counts(obj, normalized = T) %>% 
    as.data.frame(row.names = NULL) %>% 
    rownames_to_column("ensembl_gene_id") %>% 
    dplyr::filter(ensembl_gene_id %in% ensembl_goi) %>% 
    pivot_longer(cols = 2:ncol(.), values_to = "counts", names_to = "sample")
  
  md = colData(obj) %>% 
    as.data.frame(row.names = NULL) %>% 
    rownames_to_column("sample")
  
  out = left_join(genes, md, by = "sample") %>% 
    left_join(., gene_conv, by = "ensembl_gene_id") # for visualization with symbols
  
  if(df_shape == "wide") #may cause issues with overlapping gene names, ugly formatting; not recommended.
  {
    out = out %>% 
      dplyr::select(-ensembl_gene_id) %>% 
      pivot_wider(names_from = external_gene_name, values_from = counts)
  }
  
  return(out)
}
    