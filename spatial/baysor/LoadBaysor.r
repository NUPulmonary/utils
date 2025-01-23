LoadBaysor = function(data.dir, 
         cell_feature_dir, #for alternate Baysor output, normally data.dir/cell_feature_matrix/
         fov = 'fov', 
         assay = 'Xenium') {
  data <- ReadBaysor(
    data.dir = data.dir,
    cell_feature_dir =  cell_feature_dir,
    type = c("centroids", "segmentations"),
  )
  
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
  if("Blank Codeword" %in% names(data$matrix))
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Blank Codeword"]])
  else
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  
  xenium.obj[[fov]] <- coords
  return(xenium.obj)
}

ReadBaysor = function(
    data.dir,
    cell_feature_dir, #for alternate Baysor output, normally data.dir/cell_feature_matrix/
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
        matrix <- suppressWarnings(Read10X(data.dir = cell_feature_dir))
        pmtx(type = "finish")
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
          cell_info <- as.data.frame(data.table::fread(file.path(data.dir, "cells.csv.gz")))
        } else {
          cell_info <- read.csv(file.path(data.dir, "cells.csv.gz"))
        }
        cell_centroid_df <- data.frame(
          x = cell_info$x_centroid,
          y = cell_info$y_centroid,
          cell = cell_info$cell_id,
          stringsAsFactors = FALSE
        )
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
        if (has_dt) {
          cell_boundaries_df <- as.data.frame(data.table::fread(file.path(data.dir, "cell_boundaries.csv.gz")))
        } else {
          cell_boundaries_df <- read.csv(file.path(data.dir, "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
        }
        names(cell_boundaries_df) <- c("cell", "x", "y")
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
          tx_dt <- as.data.frame(data.table::fread(file.path(data.dir, "transcripts.csv.gz")))
          transcripts <- subset(tx_dt, qv >= mols.qv.threshold)
        } else {
          transcripts <- read.csv(file.path(data.dir, "transcripts.csv.gz"))
          transcripts <- subset(transcripts, qv >= mols.qv.threshold)
        }
        
        df <-
          data.frame(
            x = transcripts$x_location,
            y = transcripts$y_location,
            gene = transcripts$feature_name,
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