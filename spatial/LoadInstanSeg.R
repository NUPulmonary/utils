#' Load InstanSeg outputs into Seurat
#' 
#' @param json_path path to the GEOJSON output from InstanSeg
#' @param fov name for the fov/image object
#' @param assay name of the assay object (defaults to "CODEX")
#' @param measurement_type mean or median expression values (default: Mean)
#' @param normalize_bitdepth scale fluorescence intensity to 0-1? Avoids issues with 8-bit and 16-bit images being used in the same dataset. (default: TRUE)
#' @param perform_normalization create "data" slot with arcsinh normalization? Set to false to perform your own normalization later. (default: TRUE)
#' @param normalization_cofactor cofactor for asinh normalization. (default: 1)
#' @return a Seurat v5 object with full spatial information
#' @import Seurat sf jsonlite dplyr tibble tidyr purrr magrittr
#' @export

LoadInstanSeg = function(json_path,
                         fov = "fov",
                         assay = "CODEX",
                         normalize_bitdepth = TRUE,
                         perform_normalization = TRUE,
                         normalization_cofactor = 1,
                         measurement_type = c("Mean", "Median")){
  
  sf::sf_use_s2(FALSE) #by default assumes lat/long data
  
  measurement_type = match.arg(measurement_type)
  raw = st_read(json_path) %>% 
    dplyr::filter(!is.na(measurements)) %>% 
    #remove invalid segmentations
    dplyr::mutate(valid = st_is_valid(geometry)) %>% 
    dplyr::filter(valid == TRUE) %>% 
    dplyr::rename(cell = id)
  
  #### Create feature matrix ####
  remove_measurement_type = function(col){
    out = gsub("Cell: ", "", col)
    out = gsub(paste(":", measurement_type), "", out)
    return(make.names(out)) }
  remove_measurement_type_cytonuc = function(col){
    out = gsub("Cytoplasm: ", "Cytoplasm_", col)
    out = gsub("Nucleus: ", "Nucleus_", out)
    out = gsub(paste(":", measurement_type), "", out)
    return(make.names(out)) }
  
  features = raw %>% 
    dplyr::select(cell, measurements) %>% 
    #This is a JSON of JSONs. Force processing.
    dplyr::mutate(cell = make.names(cell),
                  measurements = map(measurements, fromJSON),
                  measurements = map(measurements, ~ map_if(.x, is.character, as.numeric))) %>% 
    #merge into single tibble
    unnest_wider(measurements, simplify = TRUE, strict = TRUE) %>% 
    dplyr::rename(nucleus_area = `Nucleus: Area µm^2`,
                  cell_area = `Cell: Area µm^2`,
                  nucleus_circularity = `Nucleus: Circularity`,
                  cell_circularity = `Cell: Circularity`,
                  nucleus_solidity = `Nucleus: Solidity`,
                  cell_solidity = `Cell: Solidity`) %>% 
    dplyr::select(cell, nucleus_area, cell_area, 
                  nucleus_circularity, cell_circularity,
                  nucleus_solidity, cell_solidity,
                  #keep only mean or median measurements for whole cell
                  matches(paste0("^Cell: .+: ", measurement_type))) %>% 
    dplyr::rename_with(.fn = remove_measurement_type) %>% 
    column_to_rownames("cell") %>% 
    t()
  
  #grab nuclear and cytoplasmic data as well for more complex analysis
  #to add: pivot longer by nuc/cyto and then get ratio
  cytonuclear = raw %>% 
    dplyr::select(cell, measurements) %>% 
    #This is a JSON of JSONs. Force processing.
    dplyr::mutate(cell = make.names(cell),
                  measurements = map(measurements, fromJSON),
                  measurements = map(measurements, ~ map_if(.x, is.character, as.numeric))) %>% 
    #merge into single tibble
    unnest_wider(measurements, simplify = TRUE, strict = TRUE) %>% 
    dplyr::select(cell,
                  matches(paste0("^Cytoplasm: .+: ", measurement_type,
                                 "|","^Nucleus: .+: ", measurement_type))) %>% 
    dplyr::rename_with(.fn = remove_measurement_type_cytonuc) %>% 
    column_to_rownames("cell") %>% 
    t()
  
  #### Import segmentation data ####
  
  #centroid calculation
  centroids = st_centroid(raw) %>% 
    st_coordinates()
  rownames(centroids) = make.names(raw$cell)
  
  #restructure segmentation data as list of matrices
  coords = st_coordinates(raw) %>% 
    as.data.frame() %>% 
    #L2 is cell number
    dplyr::mutate(L2 = as.character(L2)) %>% 
    left_join(., data.frame(L2 = rownames(raw), cell = raw$cell)) %>% 
    dplyr::select(x = X, y = Y, cell) %>% 
    dplyr::mutate(cell = make.names(cell)) %>% 
    dplyr::relocate(cell, x, y)

  geoms = CreateFOV(
    coords = list(
      "centroids" = CreateCentroids(centroids),
      "segmentation" = CreateSegmentation(coords)
      ),
    type = c("segmentation", "centroids"),
    assay = assay
    )
  
  #### scale to bit depth ####
  if(normalize_bitdepth == TRUE)
  {
    #features
    #normalize only fluorescence rows
    fluorescence_rows = which(!grepl("_area|_circularity|_solidity", rownames(features)))
    bitdepth = ifelse(max(features[fluorescence_rows, ], na.rm = TRUE) <= 255,
                      yes = 255,
                      no = 65535)
    features[fluorescence_rows, ] = features[fluorescence_rows, ] / bitdepth
    
    #cytonuclear
    #only fluorescence rows in this one
    cytonuclear = cytonuclear / bitdepth
  }
  
  #### Construct Seurat Object ####
  obj = CreateSeuratObject(counts = features, assay = assay)
  obj[[fov]] = geoms
  obj@misc$nuclear_cytoplasmic_raw = cytonuclear
  
  #### Perform Normalization ####
  if(perform_normalization == TRUE)
  {
    features_norm = features
    
    #start with morphology measures (log)
    size_rows = which(grepl("_area", rownames(features_norm)))
    features_norm[size_rows, ] = log(features_norm[size_rows, ])
    
    #fluorescence features (asinh)
    fluorescence_rows = which(!grepl("_area|_circularity|_solidity", rownames(features_norm)))
    features_norm[fluorescence_rows, ] = asinh(features_norm[fluorescence_rows, ] / normalization_cofactor)
    #format according to Seurat oddities
    rownames(features_norm) = gsub("_", "-", rownames(features_norm))
    
    #circularity and solidity measures are already on a reasonable scale (0-1)
    #and fairly normally distributed
    
    #cytonuclear
    cytonuclear_norm = asinh(cytonuclear / normalization_cofactor)
    
    #update Seurat object
    obj@assays[[assay]]$data = features_norm
    obj@misc$nuclear_cytoplasmic_norm = cytonuclear_norm
  }
  
  return(obj)
}