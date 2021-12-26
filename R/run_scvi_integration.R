#' Function to perform SCVI integration on a Seurat object
#' 
#' @param object Seurat object after variable gene selection and normalization. Can be object itself or path to an RDS file. Note: raw gene counts are used in either case.
#' @param python_path path to desired python version for reticulate setup
#' @param batch_col name of metadata column specifying batch. Defaults to orig.ident, where each sample is its own batch.
#' @param RDS_path location to save RDS package. Ignored if save_RDS is FALSE.
#' @param use_GPU toggle to true for running on GPU
#' 
#' @return a Seurat object with new embeddings in the "SCVI" slot
#' @export
#' 
run_SCVI_integration = function(object, 
                                python_path = NA,
                                batch_col = "orig.ident",
                                RDS_path = NA,
                                use_GPU = FALSE)
{
  require(Seurat)
  require(sceasy)
  
  #set up reticulate
  require(reticulate)
  if(!is.na(python_path))
  {
    use_python(python_path, required = T)
  }
  py_config()
  
  #import python packages
  sc = import('scanpy', convert = FALSE)
  scvi = import('scvi-tools', convert = FALSE)
  
  scvi$settings$progress_bar_style = 'tqdm'
  
  #if filepath is given as input, load as an RDS file
  if(class(object) == "character")
  {
    object = readRDS(object)
  }
  
  #get variable features in case SCT was used
  object = FindVariableFeatures(object, selection.method = "vst", nfeatures = 3000)
  top3k = head(VariableFeatures(object), 3000)
  object = object[top3k]
  
  #convert to annData object
  annData = convertFormat(object, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
  #filter empty cells
  
  #create model
  scvi$model$SCVI$setup_anndata(annData, batch_key = batch_col)
  model = scvi$model$SCVI(annData)
  
  # train the model
  model$train(use_gpu = useGPU)
  
  #get latent representation and place back into Seurat object
  latent = model$get_latent_representation()
  latent = as.matrix(latent)
  rownames(latent) = colnames(object)
  object[["SCVI"]] = CreateDimReducObject(embeddings = latent, key = "SCVI_", assay = DefaultAssay(object))
  
  #finally, return as necessary
  if(!is.na(RDS_path))
  {
    saveRDS(object, RDS_path)
  }
}