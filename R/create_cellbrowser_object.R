#' Function to generate all necessary files for creating a UCSC cellbrowser object from a Seurat object
#' 
#' @param object the seurat object of interest, completely finalized
#' @param ident the ident to use for cluster identification (character name). If NA (default) uses current ident
#' @param list_col_action whether to flatten, remove, or ignore list cols
#' @param output_directory directory to ouput files. Created if necessary.
#' @export

create_cellbrowser_object = function(object, ident = NA, output_directory, list_col_action = c("flatten", "remove", "ignore"))
{
  library(Seurat)
  library(SeuratDisk)
  library(tidyverse)
  
  #reset Idents as needed
  if(!is.na(ident))
  {
    object = SetIdent(object = object, ident.use = ident)
  }
  
  #create output directory as needed
  if(!dir.exists(output_directory))
  {
    dir.create(output_directory, recursive = F)
  }
  
  #create H5AD file of counts
  outname = paste0(output_directory, "/", basename(output_directory), ".h5ad")
  SaveH5Seurat(object = object, filename = outname)
  
  #output metadata
  outname = paste0(output_directory, "/", basename(output_directory), "-metadata.csv")
  md = object@meta.data
  list_cols = which(vapply(md, FUN = is.list, FUN.VALUE = T)) #get list col indices
  if(list_col_action == "flatten")
  {
    for(col in list_cols)
    {
      md = md %>% 
        rowwise() %>% 
        dplyr::mutate(across(all_of(list_cols), function(x){ paste(x, collapse = " ") }))
    }
  } else if(list_col_action == "remove")
  {
    md = md %>% 
      dplyr::select(-all_of(list_cols))
  } else if(list_col_action == "ignore")
  {
  } else
  {
    stop("Error: incorrect list_col_action selection")
  }
    write.csv(md, outname)
    
    #output markers
    outname = paste0(output_directory, "/", basename(output_directory), "-markers.csv")
    markers = FindAllMarkers(object = object, only.pos = TRUE)
    write.csv(markers, outname)
}  