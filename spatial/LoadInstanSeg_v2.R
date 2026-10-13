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
#' @import Seurat sf jsonlite dplyr tibble tidyr purrr magrittr janitor
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
  
  raw = sf::st_read(
    json_path,
    options = "FLATTEN_NESTED_ATTRIBUTES=YES",
    quiet = TRUE
  ) %>% 
    clean_names(replace = c("µ" = "u", "iNOS" = "iNOS", "NeuN" = "NeuN",
                            "Ki67" = "Ki67", "CD3e" = "CD3e"), case = "none") %>% 
    dplyr::filter(objectType == "cell") %>%   #get rid of annotations
    dplyr::rename_with(~ gsub("measurements_", "", .x))
  
  #keep only columns with desired measurement type (e.g. Mean) and morphological measures
  keep_cols = colnames(raw)[grepl(measurement_type, colnames(raw))]
  keep_cols = c("id", "geometry", "Nucleus_Area_um_2","Nucleus_Circularity", 
                "Nucleus_Solidity", "Cell_Area_um_2","Cell_Circularity", 
                "Cell_Solidity", keep_cols)
  raw = raw[, keep_cols] %>%
    dplyr::mutate(
      across(3:ncol(.), as.numeric)
    ) %>%
    dplyr::filter(
      if_all(3:ncol(.), is.finite),
      Nucleus_Area_um_2 > 0,
      Cell_Area_um_2 > 0
    ) %>%
    dplyr::rename(cell = id)
  
  ####### Feature matrix #######
  #cell features only
  keep_cols_features = colnames(raw)[grepl(paste0("Cell_.+", measurement_type), colnames(raw))]
  keep_cols_features = c("cell", "Nucleus_Area_um_2","Nucleus_Circularity", 
                         "Nucleus_Solidity", "Cell_Area_um_2","Cell_Circularity", 
                         "Cell_Solidity", keep_cols_features)
  
  
  features = raw[, keep_cols_features] %>% 
    column_to_rownames("cell") %>% 
    dplyr::rename_with(~ str_replace(.x, paste0("_", measurement_type), "")) %>% 
    dplyr::rename_with(~ str_replace(.x, "_um_2", "")) %>% 
    dplyr::rename_with(~ str_replace(.x, "^Cell_", "")) %>% 
    #hard to avoid this
    dplyr::rename(Cell_Area = Area,
                  Cell_Circularity = Circularity,
                  Cell_Solidity = Solidity) %>% 
    dplyr::select(-geometry) %>% 
    t() %>% 
    as.matrix()
  
  ####### Nuclear/Cytoplasmic values #######
  fluorescence_rows = which(!grepl("_Area|_Circularity|_Solidity", rownames(features)))
  bitdepth = ifelse(max(features[fluorescence_rows, ], na.rm = TRUE) <= 255,
                    yes = 255,
                    no = 65535)
  pseudocount_raw = bitdepth / 65535
  
  keep_cols_cytonuclear = colnames(raw)[grepl(paste0("(Nucleus|Cytoplasm)_.+", measurement_type), colnames(raw))]
  keep_cols_cytonuclear = c("cell", keep_cols_cytonuclear)
  cytonuclear = raw[, keep_cols_cytonuclear] %>% 
    dplyr::rename_with(~ str_replace(.x, paste0("_", measurement_type), "")) %>% 
    as.data.frame() %>% 
    dplyr::select(-geometry) %>% 
    #infinite/NA values = no nucleus or cytoplasm
    dplyr::filter(if_all(where(is.numeric), is.finite)) %>% 
    #now create ratios
    pivot_longer(
      cols = matches("^(Nucleus|Cytoplasm)"),
      names_to = c(".value", "protein"),
      names_pattern = paste0("^(Nucleus|Cytoplasm)_(.+)")
    ) %>% 
    dplyr::mutate(Nucleus = as.numeric(Nucleus),
                  Cytoplasm = as.numeric(Cytoplasm),
                  #pseudocount of a sort to avoid Inf
                  NC_ratio = (Nucleus + pseudocount_raw) / (Cytoplasm + pseudocount_raw))

  #### Import segmentation data ####
  #centroid calculation
  centroids = st_centroid(raw) %>% 
    st_coordinates()
  rownames(centroids) = raw$cell
  centroids_df = data.frame(
    cell = raw$cell,
    x = centroids[,1],
    y = centroids[,2]
  )
  
  #segmentation coordinates
  coords_list = lapply(seq_len(nrow(raw)), function(i) {
    
    cur_geom = raw$geometry[i]
    cur_cell = raw$cell[i]
    
    # extract coordinates for this geometry only
    cur_coords = st_coordinates(cur_geom) %>%
      as.data.frame()
    
    # retain only x/y
    cur_coords = cur_coords[, c("X", "Y")]
    colnames(cur_coords) = c("x", "y")
    
    # attach cell identity explicitly
    cur_coords$cell = cur_cell
    
    # order columns
    cur_coords = cur_coords[, c("cell", "x", "y")]
    
    return(cur_coords)
  })
  
  coords = dplyr::bind_rows(coords_list)
  
  geoms = CreateFOV(
    coords = list(
      "centroids" = CreateCentroids(centroids_df, nsides = Inf, radius = 70),
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
    pseudocount = 1 / 65535
    
    features[fluorescence_rows, ] = features[fluorescence_rows, ] / bitdepth
    
    #cytonuclear
    #only fluorescence rows in this one
    cytonuclear = cytonuclear %>% 
      dplyr::mutate(Nucleus = Nucleus / bitdepth,
                    Cytoplasm = Cytoplasm / bitdepth,
                    #pseudocount of a sort to avoid Inf
                    NC_ratio = (Nucleus + pseudocount) / (Cytoplasm + pseudocount))
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
    size_rows = which(grepl("_Area", rownames(features_norm)))
    features_norm[size_rows, ] = log(features_norm[size_rows, ])
    
    #fluorescence features (asinh)
    fluorescence_rows = which(!grepl("_Area|_Circularity|_Solidity", rownames(features_norm)))
    features_norm[fluorescence_rows, ] = asinh(features_norm[fluorescence_rows, ] / normalization_cofactor)
    #format according to Seurat oddities
    rownames(features_norm) = gsub("_", "-", rownames(features_norm))
    
    #circularity and solidity measures are already on a reasonable scale (0-1)
    #and fairly normally distributed
    
    #cytonuclear
    cytonuclear_norm = cytonuclear %>% 
      dplyr::mutate(Nucleus = asinh(Nucleus / normalization_cofactor),
                    Cytoplasm = asinh(Cytoplasm / normalization_cofactor))
    
    #update Seurat object
    obj@assays[[assay]]$data = features_norm
    obj@misc$nuclear_cytoplasmic_norm = cytonuclear_norm
  }
  return(obj)
}
