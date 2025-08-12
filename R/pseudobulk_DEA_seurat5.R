#' Function to separate single-cell data by a factor of interest and create  pseudo-bulk datasets by combining matrices for each factor followed by DEA with DESeq2   
#' Note: currently supports only single-factor designs
#' Updated 2025-08-12 to support Seurat 5+ RAG
#' 
#' @param object seurat object
#' @param assay the assay from the object to pseudobulk. Defaults to RNA.
#' @param sample_metadata metadata mapping samples to factors. Rownames are sample names.
#' @param design design object (NOT STRING) for DESeq
#' @param skip if splitting cells, which factors should be skipped? Defaults to none.
#' @param sample_factor name of the sample column (i.e. group.by). Defaults to orig.ident
#' @param cell_factor name of the cell-type column
#' @param organism organism in ensembl format, e.g. mmusculus
#' @param gene_var type of gene ID used, in biomart format e.g. ensembl_gene_id
#' @param out_dir directory to output results
#' @param out_prefix file prefix for results files
#' @param sort_direction direction to sort comparisons: alphabetical "ascending" or reverse "descending"
#' @param min_cells whether or not to override the minimum number of cells/samples
#' @param cores cores to run in parallel
#' @param fit_type to control dispersion fitting
#' @param genome_prefix regular expression of genome prefixe(s) on gene names for removal and better binding downstream
#' @import Seurat DESeq2 biomaRt future BiocParallel tibble
#' @return A list of DESeq2 dge objects. All relevant CSVs, PDFs of plots, and RDS files are saved to the directory specified.
#' @export
  
pseudobulk_DEA = function(object, 
                   assay = "RNA",
                   design, 
                   skip = NA,
                   sample_factor = "orig.ident",
                   cell_factor, 
                   organism = "mmusculus", 
                   gene_var,
                   out_dir, 
                   out_prefix, 
                   sort_direction = c("factor_order", "descending", "ascending"), 
                   min_cells = 50, 
                   cores = 1,
                   fit_type = "parametric",
                   genome_prefix = NA, 
                   cell_mappings = NULL) 
{
  
  if(!("Seurat" %in% .packages()))
  {
    library(Seurat)
  }
  if(!("DESeq2" %in% .packages()))
  {
    library(DESeq2)
  }
  if(!("biomaRt" %in% .packages()))
  {
    library(biomaRt)
  }
  if(!("future" %in% .packages()))
  {
    library(future)
  }
  if(!("BiocParallel" %in% .packages()))
  {
    library(BiocParallel)
  }
  if(!("tibble" %in% .packages()))
  {
    library(tibble)
  }
  
  #calling formula in function causes massive problems. Avoid by doing externally!
  if(class(design) != "formula")
  {
    stop("design paramater must be a formula object. Calling formula in function causes massive problems.")
  }
  
  #extract metadata
  comp_col = as.character(design)[2]
  sample_metadata = unique(object@meta.data[, c(sample_factor, comp_col)])
  
   register(BPPARAM = MulticoreParam(cores))
   if(!dir.exists(out_dir))
      dir.create(out_dir)
   setwd(out_dir)
   
   sort_direction = match.arg(sort_direction)
   
   #make mart for conversion
   mart = useMart("ensembl", paste0(organism, "_gene_ensembl"))
   conv = getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name"), mart = mart)
   
   #change metadata to match seurat standards
   object@meta.data[, sample_factor] = gsub("_", "-", object@meta.data[, sample_factor])
   sample_metadata[, sample_factor] = gsub("_", "-", sample_metadata[, sample_factor])
   
   #cells per sample per celltype
   sample_data = table(object@meta.data[, sample_factor],
                       object@meta.data[, cell_factor])
   
   #samples per celltype (in case of missing cells for some samples)
   n_data = colSums(sample_data > min_cells)
   
   #get complete summarized counts by sample, celltype
   all_pb_counts = AggregateExpression(object,
                                       assays = assay,
                                       group.by = c(cell_factor, sample_factor), #this order allows for numeric sample names
                                       return.seurat = FALSE) #this forces no normalization. Essential for DESeq2.
   all_pb_counts = as.matrix(all_pb_counts[[assay]])
   
   #remove celltype x sample combos with insufficient cells
   combo_counts = table(paste(object@meta.data[, cell_factor], 
                              object@meta.data[, sample_factor],
                              sep = "_"))
   good_combos = names(combo_counts)[combo_counts > min_cells]
   all_pb_counts = all_pb_counts[, good_combos]
   
   # determine all possible pairwise comparisons (and avoid repeating this)
   all_comps = combn(levels(object@meta.data[, comp_col]), 2, simplify = F)
   all_comps = lapply(all_comps, function(comp){
     if(sort_direction == "factor_order")
     {
       comp = factor(comp, levels = levels(object@meta.data[, comp_col]))
       comp = sort(comp)
     }
     else if(sort_direction == "ascending")
     {
       comp = sort(comp)
     } else if(sort_direction == "descending")
     {
       comp = sort(comp, decreasing = T)
     }
     comp = as.character(comp) #get in format for results()
     return(comp) })
   
   celltypes = unique(object@meta.data[ , cell_factor])
   celltypes = setdiff(celltypes, skip)
   dges = lapply(celltypes, function(cell){
     
     message(cell)
     
     cur_mat = all_pb_counts[, grepl(cell, colnames(all_pb_counts))]
     cur_md = sample_metadata
     rownames(cur_md) = paste(cell, cur_md[, sample_factor], sep = "_") #in case sample is numeric, this is easier
     cur_md = cur_md[colnames(cur_mat), ]
     
     #handle cases with insufficient n
     des = tryCatch({ DESeqDataSetFromMatrix(countData = cur_mat, colData = cur_md, design = design) },
                    error = function(msg){
                      warning("Skipped due to insufficient n")
                      return(NULL) })
     if(is.null(des))
     {
       return(NULL)
     }
     dge = DESeq(des, parallel = T, fitType = fit_type)
     
     #output DGE. Safe for symbols.
     out_path = paste(out_prefix, gsub("\\W", replacement = "-", cell), "subset", "des.rds", sep = "_")
     saveRDS(dge, out_path)
     
     
      for(comp in all_comps)
      {
         res = as.data.frame(results(dge, contrast = c(comp_col, comp[1], comp[2]), alpha = 0.05, parallel = TRUE))
         res = rownames_to_column(res, var = gene_var)
         if(!is.na(genome_prefix) && gene_var == "external_gene_name")
         {
           res[, gene_var] = gsub(genome_prefix, "", res[, gene_var])
         }
         if(gene_var == "ensembl_gene_id")
         {
            res$ensembl_gene_id = ifelse(grepl("\\.", res$ensembl_gene_id),
                                         yes = substr(res$ensembl_gene_id, 1, (regexpr("\\.", res$ensembl_gene_id) - 1)),
                                         no = res$ensembl_gene_id) #remove version numbers
         }
         res = merge(res, conv, all.x = TRUE, all.y = FALSE)
         
         comp_name = paste(comp[1], "vs", comp[2], sep = "_")
         out_path = paste(out_prefix, gsub("\\W", replacement = "-", cell), "subset", comp_name, "dge.csv", sep = "_")
         write.csv(res, out_path, row.names = FALSE)
      }
     
     return(dge)
     
   })
   
   names(dges) = celltypes
   return(dges)
   
}
