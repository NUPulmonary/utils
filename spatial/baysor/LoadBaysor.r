LoadBaysor = function(baysor_dir, #for alternate Baysor output, normally data.dir/cell_feature_matrix/
         remove_bad_codewords = TRUE, #remove unassigned or deprecated
         fov = 'fov', 
         assay = 'Xenium') {
  data <- ReadBaysor(
    baysor_dir =  baysor_dir,
    type = c("centroids", "segmentations"),
  )
  
  #baysor has some very small discrepancies in cells between outputs (like 2 total cells)
  #subset down to just intersection and get into same order
  common_cells = intersect(intersect(colnames(data$matrix[["Gene Expression"]]), 
                                     data$centroids$cell),
                           data$segmentations$cell)
  data$matrix[["Gene Expression"]] = data$matrix[["Gene Expression"]][, common_cells]
  data$centroids = data$centroids[common_cells, ]
  data$segmentations = subset(data$segmentations, cell %in% common_cells)
  data$segmentations$cell = factor(data$segmentations$cell, levels = rownames(data$centroids))
  data$segmentations = data$segmentations[order(data$segmentations$cell), , drop = FALSE]
  
  segmentations.data <- list(
    "centroids" = CreateCentroids(data$centroids),
    "segmentation" = CreateSegmentation(data$segmentations)
  )
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = data$microns,
    assay = assay
  )
  
  xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
  
  #all lost with Baysor AFAIK
  # if("Blank Codeword" %in% names(data$matrix))
  #   xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Blank Codeword"]])
  # else
  #   xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
  # xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  # xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  
  xenium.obj[[fov]] <- coords
  return(xenium.obj)
}

ReadBaysor = function(
    baysor_dir, #for alternate Baysor output, normally data.dir/cell_feature_matrix/
    remove_bad_codewords = TRUE, #remove unassigned or deprecated
    outs = c("matrix", "microns"),
    type = "centroids",
    mols.qv.threshold = 20
) {
  # Argument checking
  type <- match.arg(
    arg = type,
    choices = c("centroids", "segmentations"),
    several.ok = TRUE
  )
  
  outs <- match.arg(
    arg = outs,
    choices = c("matrix", "microns"),
    several.ok = TRUE
  )
  
  outs <- c(outs, type)
  
  has_dt <- requireNamespace("data.table", quietly = TRUE) && requireNamespace("R.utils", quietly = TRUE)
  
  data <- sapply(outs, function(otype) {
    switch(
      EXPR = otype,
      'matrix' = {
        pmtx <- progressr::progressor()
        pmtx(message = 'Reading counts matrix', class = 'sticky', amount = 0)
        #note!! Baysor currently does not split the matrix into categories, but may at some point. 
        #Fix next 8 lines accordingly.
        matrix <- suppressWarnings(Read10X(data.dir = paste0(baysor_dir, "/baysor_mtx")))
        if(remove_bad_codewords == TRUE)
        {
          good_codewords = rownames(matrix)[!(grepl("UnassignedCodeword|DeprecatedCodeword|BlankCodeword|Negative Control Codeword|Negative Control Probe", rownames(matrix)))]
          matrix = matrix[good_codewords, ]
        }
        pmtx(type = "finish")
        matrix <- list("Gene Expression" = matrix) #here as well
        matrix
      },
      'centroids' = {
        pcents <- progressr::progressor()
        pcents(
          message = 'Loading cell centroids',
          class = 'sticky',
          amount = 0
        )
        if (has_dt) {
          cell_info <- as.data.frame(data.table::fread(paste0(baysor_dir, "/segmentation_cell_stats.csv")))
        } else {
          cell_info <- read.csv(paste0(baysor_dir, "/segmentation_cell_stats.csv"))
        }
        cell_centroid_df <- data.frame(
          #note that all 3 are different in segmentation_cell_stats.csv vs cells.csv.gz
          x = cell_info$x,
          y = cell_info$y,
          #Baysor barcodes.tsv.gz does this for some reason. Need to do this to successfully link.
          cell =  paste0("cell_", cell_info$cell),
          stringsAsFactors = FALSE
        )
        rownames(cell_centroid_df) = cell_centroid_df$cell
        pcents(type = 'finish')
        cell_centroid_df
      },
      'segmentations' = {
        psegs <- progressr::progressor()
        psegs(
          message = 'Loading cell segmentations',
          class = 'sticky',
          amount = 0
        )
        
        # load cell boundaries
        cell_boundaries_df = parse_Baysor_JSON(paste0(baysor_dir, "/segmentation_polygons_2d.json"))
        psegs(type = "finish")
        cell_boundaries_df
      },
      'microns' = {
        pmicrons <- progressr::progressor()
        pmicrons(
          message = "Loading molecule coordinates",
          class = 'sticky',
          amount = 0
        )
        
        # molecules
        if (has_dt) {
          tx_dt <- as.data.frame(data.table::fread(paste0(baysor_dir, "/segmentation.csv")))
          transcripts <- subset(tx_dt, qv >= mols.qv.threshold) #qv keeps same name
        } else {
          transcripts <- read.csv(paste0(baysor_dir, "/segmentation.csv"))
          transcripts <- subset(transcripts, qv >= mols.qv.threshold)
        }
        
        df <-
          data.frame(
            #note name changes from segmentation.csv vs transcripts.csv.gz
            x = transcripts$x,
            y = transcripts$y,
            gene = transcripts$gene,
            stringsAsFactors = FALSE
          )
        pmicrons(type = 'finish')
        df
      },
      stop("Unknown Xenium input type: ", otype)
    )
  }, USE.NAMES = TRUE)
  return(data)
}

parse_Baysor_JSON = function(filepath){
  json = jsonlite::fromJSON(filepath, flatten = TRUE)
  json = json$features
  #extract lists of x, y coords
  json$x = lapply(json$geometry.coordinates, function(row){
    return(row[,,1]) })
  json$y = lapply(json$geometry.coordinates, function(row){
    return(row[,,2]) })
  
  #now pivot to single entries (long form)
  json = tidyr::unnest_longer(json, col = c(x, y))
  
  #simplify output
  out = data.frame(cell = json$id, x = json$x, y = json$y)
  #Baysor barcodes.tsv.gz does this for some reason. Need to do this to successfully link.
  out$cell = paste0("cell_", out$cell)
  return(out)
}
