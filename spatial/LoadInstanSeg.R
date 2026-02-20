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
  centroids = st_centroid(raw$geometry) %>% 
    st_coordinates()
  rownames(centroids) = make.names(raw$cell)
  
  #restructure segmentation data as list of matrices
  coords = st_coordinates(raw$geometry) %>% 
    as.data.frame() %>% 
    #L3 is cell number
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
  
  #### Construct Seurat Object ####
  obj = CreateSeuratObject(counts = features, assay = assay)
  obj[[fov]] = geoms
  obj@misc$nuclear_cytoplasmic = cytonuclear
  return(obj)
}
