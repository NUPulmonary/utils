#' Function to perform SCVI integration on a Seurat object
#' 
#' @param object Seurat object after variable gene selection and normalization. Can be object itself or path to an RDS file. Note: raw gene counts are used in either case.
#' @param python_path path to desired python version for reticulate setup
#' @param random_seed random seed for SCVI
#' @param project_path path to project directory for loading R environment. If NA, does not use renv.
#' @param batch_col name of metadata column specifying batch. Defaults to orig.ident, where each sample is its own batch.
#' @param RDS_path location to save RDS package. Ignored if save_RDS is FALSE.
#' @param use_GPU toggle to true for running on GPU
#' @param hvgs number of highly variable genes to use in the model. Default: 3000.
#' @param n_layers Number of hidden layers used for encoder and decoder NNs
#' @param dropout_rate Dropout rate for neural networks
#' 
#' @return a Seurat object with new embeddings in the "SCVI" slot
#' @export
#' 
run_SCVI_integration = function(object, 
                                python_path = NA,
                                project_path = NA,
                                random_seed = 12345,
                                batch_col = "orig.ident",
                                n_epochs = NA,
                                RDS_path = NA,
                                use_GPU = FALSE,
                                hvgs = 3000,
                                n_layers = 1,
                                dropout_rate = 0.1,
                                early_stopping = FALSE,
                                n_latent = 10)
{
  #if required, set up R environment
  if(!is.na(project_path))
  {
    #need renv installed outside of project folder for this to work
    if(!("renv" %in% installed.packages()))
    {
      install.packages("renv")
    }
    require(renv)
    renv::load(project_path)
  }
  
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
  scvi = import('scvi', convert = FALSE)
  
  scvi$settings$progress_bar_style = 'tqdm'
  scvi$settings$seed = as.integer(random_seed) #causes scvi error unless cast
  
  #if filepath is given as input, load as an RDS file
  if(class(object) == "character")
  {
    object = readRDS(object)
  }
  
  #get variable features (assume already calculated with method of choice)
  top_n = head(VariableFeatures(object), hvgs)
  #subsetting in Seurat 5 is extremely weird. Handle as necessary
  if(substring(packageVersion("Seurat"), 1, 1) <= 4) #version 4 and below
  {
    object_scvi = object[top_n, ]
  } else #version 5 and up
  {
    object_scvi = object
    object_scvi[[DefaultAssay(object_scvi)]] = subset(object_scvi[[DefaultAssay(object_scvi)]], features = top_n)
    if(DefaultAssay(object_scvi) != "RNA")
    {
      object_scvi[["RNA"]] = subset(object_scvi[["RNA"]], features = top_n)
    }
  }
  
  #convert to annData object
  annData = convertFormat(object_scvi, 
                          from="seurat", 
                          to="anndata", 
                          main_layer="counts", 
                          drop_single_values=FALSE)
  
  #create model
  scvi$model$SCVI$setup_anndata(annData, batch_key = batch_col)
  model = scvi$model$SCVI(annData, dropout_rate = as.double(dropout_rate), n_layers = as.integer(n_layers))
  
  #if supplied, set number of epochs
  
  
  # train the model
  if(is.na(n_epochs))
  {
    model$train(use_gpu = use_GPU, early_stopping = as.logical(early_stopping))
  } else
  {
    model$train(use_gpu = use_GPU, max_epochs = as.integer(n_epochs), early_stopping = as.logical(early_stopping))
  }
  
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