#' Process a Seurat5 object for cellbrowser upload
#' 
#' @param obj a Seurat object (v5) or path to a Seurat v5 RDS file
#' @param fov which FOV to export for image "reduction". If there are multiple FOVs, recommend running concatenate_fovs() first. If NULL (default), no image reduction is included.
#' @param assay assay to pull counts from. Defaults to RNA.
#' @param metadata_cols names of metadata columns to include. If NULL (default) all columns are included.
#' @param dataset_name name of the dataset for cellbrowser page
#' @param short_name short, URL-safe name for dataset (shown in url on cellbrowser)
#' @param cluster_col name of column to use for default clustering. Defaults to seurat_clusters.
#' @param outdir directory for export of processed files
#' @param compress_directory whether or not to tar the entire directory after completion (for upload)
#' @param delete_files whether or not to delete constituent files after compression (ignored if compress_directory is FALSE)
#' @return Nothing; just exports files
#' @import Seurat readr tibble R.utils
#' @export

cellbrowser_export_seurat5 = function(obj,
                                      fov = NULL,
                                      assay = "RNA",
                                      metadata_cols = NULL,
                                      dataset_name,
                                      short_name,
                                      cluster_col = "seurat_clusters",
                                      outdir,
                                      compress_directory = TRUE,
                                      delete_files = FALSE){
  
  if(!("Seurat" %in% .packages()))
  {
    library(Seurat)
  }
  #need readr for faster output (matters with files this large)
  if(!("readr" %in% .packages()))
  {
    library(readr)
  }
  if(!("tibble" %in% .packages()))
  {
    library(tibble)
  }
  if(!("R.utils" %in% .packages()))
  {
    library(R.utils)
  }
  
  if(!dir.exists(outdir))
  {
    dir.create(outdir, recursive = TRUE)
  }
  setwd(outdir)
  
  #start making config file
  write_lines(x = c(paste0("name = \"", short_name, "\""),
                    "tags = [\"10x\"]",
                    paste0("shortLabel = \"", dataset_name, "\""),
                    "exprMatrix = \"exprMatrix.tsv.gz\"",
                    "geneIdType = \"symbol\"",
                    "meta = \"meta.tsv\"",
                    "markers = \" markers.tsv\"",
                    #avoid misuse of factor cols that look numeric
                    paste0("enumFields = [\"", paste(colnames(obj@meta.data)[unlist(lapply(obj@meta.data, is.factor), use.names = FALSE)], collapse = "\", \""), "\"]"),
                    paste0("clusterField = \"", cluster_col, "\""),
                    paste0("labelField = \"", cluster_col, "\"")),
              file = file("cellbrowser.conf"), 
              append = FALSE)
  
  #load RDS, as necessary
  if(class(object) == "character")
  {
    obj = readRDS(obj)
  }

  # UMAP
  message("Writing UMAP reduction data")
  umap_out = rownames_to_column(as.data.frame(obj@reductions$umap@cell.embeddings), "cellName")
  colnames(umap_out) = c("cellName", "x", "y")
  write_tsv(umap_out, "umap_coords.tsv")
  
  #Expression matrix
  message("Writing expression data")
  out_mat = as.data.frame(obj@assays[[assay]]@layers$counts, row.names = rownames(obj))
  colnames(out_mat) = colnames(obj)
  write_tsv(out_mat, "exprMatrix.tsv")
  gzip("exprMatrix.tsv")
  
  # Metadata
  message("Writing metadata")
  md_out = rownames_to_column(obj@meta.data, "cellName")
  if(!is.null(metadata_cols))
  {
    md_out = md_out[, metadata_cols]
  }
  write_tsv(md_out, "meta.tsv")
  
  # Markers
  message("Writing marker data")
  if(is.null(obj@misc$markers))
  {
    obj@misc$markers = FindAllMarkers(obj)
  }
  write_tsv(obj@misc$markers[, c("cluster", "gene", "p_val_adj")], "markers.tsv")
  
  # Spatial coordinates (as necessary)
  if(!(is.null(fov)))
  {
    message("Writing spatial data")
    #make coordinates file
    fov_out = obj@images[[fov]]$centroids@coords
    rownames(fov_out) = obj@images[[fov]]$centroids@cells
    fov_out = rownames_to_column(as.data.frame(fov_out), "cellName")
    write_tsv(fov_out, "spatial_coords.tsv")
    
    #update config
    write_lines(x = "coords=[
    {
            \"file\":\"umap_coords.tsv\", 
            \"flipY\" : True, # R files need to be flipped on the Y-axis
            \"shortLabel\":\"UMAP\"
    },
    {
            \"file\":\"spatial_coords.tsv\", 
            \"flipY\" : True, # R files need to be flipped on the Y-axis
            \"shortLabel\":\"Spatial\", 
    },
    ]",
                file = file("cellbrowser.conf"), 
                append = TRUE)
  } else
  {
    #if not, include only UMAP
    write_lines(x = "coords=[
    {
            \"file\":\"umap_coords.tsv\", 
            \"flipY\" : True, # R files need to be flipped on the Y-axis
            \"shortLabel\":\"UMAP\"
    }
    ]",
                file = file("cellbrowser.conf"), 
                append = TRUE)
  }
  
  if(compress_directory == TRUE)
  {
    message("Compressing directory")
    tar(tarfile = paste0(short_name, ".tar.gz"), 
        files = ".", 
        compression = "gzip")
    
    if(delete_files == TRUE)
    {
      message("Cleaning up")
      file.remove("umap_coords.tsv")
      file.remove("exprMatrix.tsv.gz")
      file.remove("meta.tsv")
      file.remove("markers.tsv")
      file.remove("cellbrowser.conf")
      if(!(is.null(fov)))
      {
        file.remove("spatial_coords.tsv")
      }
    }
  }
}