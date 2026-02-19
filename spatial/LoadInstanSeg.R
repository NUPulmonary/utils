#' Load InstanSeg outputs into Seurat
#' 
#' @param json_path path to the GEOJSON output from InstanSeg
#' @param fov name for the fov/image object
#' @param assay name of the assay object (defaults to "CODEX")
#' @param measurement_type mean or median expression values (default: Mean)
#' @return a Seurat v5 object with full spatial information
#' @import Seurat sf jsonlite dplyr tibble tidyr purrr magrittr s2
#' @export

LoadInstanSeg = function(json_path,
                         fov = "fov",
                         assay = "CODEX",
                         measurement_type = c("Mean", "Median")){
  
  measurement_type = match.arg(measurement_type)
  raw = st_read(json_path) %>% 
    dplyr::filter(!is.na(measurements)) %>% 
    #remove invalid segmentations
    dplyr::mutate(valid = st_is_valid(geometry)) %>% 
    dplyr::filter(valid == TRUE) %>% 
    dplyr::rename(cell = id)
  
  #### Create feature matrix ####
  
  remove_measurement_type = function(col){
    return(make.names(gsub(paste(":", measurement_type), "", col)))}
  
  features = raw %>% 
    dplyr::select(cell, measurements) %>% 
    #This is a JSON of JSONs. Force processing.
    dplyr::mutate(cell = make.names(cell),
                  measurements = map(measurements, fromJSON),
                  measurements = map(measurements, ~ map_if(.x, is.character, as.numeric))) %>% 
    #merge into single tibble
    unnest_wider(measurements, simplify = TRUE, strict = TRUE) %>% 
    dplyr::rename(area = `Area µm^2`,
                 circularity = Circularity,
                 solidity = Solidity) %>% 
    dplyr::select(cell, area, circularity, solidity, 
                  #keep only mean or median measurements
                  ends_with(measurement_type)) %>% 
    dplyr::rename_with(.fn = remove_measurement_type) %>% 
    column_to_rownames("cell") %>% 
    t()
  
  #### Import segmentation data ####
  
  #centroid calculation
  centroids = st_centroid(raw$geometry) %>% 
    st_coordinates()
  rownames(centroids) = make.names(raw$cell)
  
  #restructure segmentation data as list of matrices
  coords = st_coordinates(raw$geometry) %>% 
    as.data.frame() %>% 
    #L3 is cell number
    dplyr::mutate(L3 = as.character(L3)) %>% 
    left_join(., data.frame(L3 = rownames(raw), cell = raw$cell)) %>% 
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
  
  #### Construct Seurat Object ####
  obj = CreateSeuratObject(counts = features, assay = assay)
  obj[[fov]] = geoms
  return(obj)
}