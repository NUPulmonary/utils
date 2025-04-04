library(arrow)      
library(dplyr)       
library(tidyr)      
library(Seurat)
library(magrittr)
##################################################
## Purpose: Takes Xenium ROI Output directory and performs nucleus filtering to form a Seurat Object
##          Filtered to Nucleus Expression. Nucleus Filtered RDS will be saved directly to the ROI output directory
##
## Variables:
##      ROI_dir: Output directory bundle for the ROI/sample you wish to nucleus filter, do not add a training "/" to the path
##      qv: Transcript qv score used to filter out all transcripts falling below such cutoff, default 10x cutoff is 20
##
##################################################
create_and_save_nucleus_filtering <- function(ROI_dir, qv) {
  transcripts <- read_parquet(paste0(ROI_dir,"/transcripts.parquet"))
  message("Transcript Matrix Loaded")
  all_transcripts <- transcripts[transcripts$qv > qv, ]
  
  # Find transcripts that overlap a nucleus
  nuc_transcripts <- all_transcripts[all_transcripts$overlaps_nucleus == "1", ]
  
  # Create cell x gene dataframe
  nuc_transcripts <- as.data.frame(table(nuc_transcripts$cell_id, 
                                         nuc_transcripts$feature_name))
  names(nuc_transcripts) <- c("cell_id", "feature_name", "Count")
  nuc_transcripts <- nuc_transcripts %>% pivot_wider(names_from = "feature_name", values_from = "Count")
  message("Nucleus Transcript Matrix Created")
  # Get blanks count per nucleus
  blank_nuc_ids <- nuc_transcripts$cell_id
  blank_nuc_mat <- nuc_transcripts[, grep("BLANK", 
                                          colnames(nuc_transcripts))]
  blank_nuc_counts <- as.data.frame(rowSums(blank_nuc_mat))
  blank_nuc_counts$cell_id <- blank_nuc_ids
  
  # Remove negative controls and convert to cell x gene matrix
  nuc_transcripts <- nuc_transcripts[, grep("NegControl", 
                                            colnames(nuc_transcripts), 
                                            invert = TRUE)]
  nuc_transcripts <- nuc_transcripts[, grep("BLANK", 
                                            colnames(nuc_transcripts), 
                                            invert = TRUE)]
  nuc_transcripts <- nuc_transcripts[, grep("Unassigned", 
                                            colnames(nuc_transcripts), 
                                            invert = TRUE)]
  keep_cells <- nuc_transcripts$cell_id
  nuc_transcripts <- as.data.frame(nuc_transcripts)
  rownames(nuc_transcripts) <- keep_cells
  nuc_transcripts <- nuc_transcripts[, -1]
  nuc_transcripts <- as.matrix(t(nuc_transcripts))
  message("Nucleus Transcript Matrix Filtered")
  options(Seurat.object.assay.version = "v3")
  nucleus_object <- CreateSeuratObject(counts = nuc_transcripts,assay="RNA", project = "Nucleus Filtering")
  message("Nucleus Transcript Seurat Object Created, Saving Now")
  saveRDS(nucleus_object, file = paste0(ROI_dir,"/nucleus_counts.rds"))
}