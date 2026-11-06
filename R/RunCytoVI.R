# Fit CytoVI to Phenocycler measurements and add a Seurat reduction.
# Source this file in R, then call run_cytovi_rds(input_rds, output_rds, batch_col).
# Requires SeuratObject >= 5, reticulate, and a Python environment containing
# scvi-tools (with scvi.external.CYTOVI), anndata, numpy, and pandas.
# Input: an RDS with a CODEX assay and sample metadata. Output: a new RDS with
# reduction "CytoVI"; the original assays, metadata, identities, and FOVs persist.
# Measurements must already be transformed: reuse scale.data, or pool-scale data.
# All available features, including DAPI and morphology, are eligible by default.

#' Prepare a cell-by-feature matrix and aligned sample annotations for CytoVI.
#'
#' @param expression Feature-by-cell matrix from one Seurat expression layer.
#' @param metadata Cell metadata data.frame, with cell names as row names.
#' @param cells Unique cell names, in the desired output order.
#' @param batch_col Metadata column identifying the samples to integrate.
#' @param features Optional feature names; NULL uses all features in the layer.
#' @param exclude_features Feature names to omit; defaults to character() (none).
#' @param scale_input Whether to z-score each feature across all cells together.
#' @param exclude_morphology Whether to omit the six Cell/Nucleus morphology
#'   measurements created by LoadInstanSeg; defaults to FALSE.
#' @return List containing x (cells by features), obs (aligned batch metadata),
#'   feature_ranges (per-feature minimum/maximum), and constant_features
#'   (features dropped for having zero variance).
#' @details Morphology matching accepts hyphens or underscores in feature names.
#'   No feature categories are excluded unless requested by the caller.
#'   Never applies arcsinh again. Stops on missing cells, features, batch labels,
#'   nonfinite measurements, or insufficient variation. Does not modify inputs.
prepare_cytovi_input <- function(expression, metadata, cells, batch_col,
                                 features = NULL, exclude_features = character(),
                                 scale_input = FALSE, exclude_morphology = FALSE) {
  if (length(exclude_morphology) != 1L || !is.logical(exclude_morphology) ||
      is.na(exclude_morphology)) {
    stop("exclude_morphology must be TRUE or FALSE.")
  }
  if (length(cells) < 3L || anyDuplicated(cells) || anyNA(cells)) {
    stop("Provide at least three cells with unique, nonmissing cell names.")
  }
  if (anyDuplicated(colnames(expression)) ||
      anyDuplicated(rownames(expression)) ||
      !all(cells %in% colnames(expression)) ||
      !all(cells %in% rownames(metadata))) {
    stop("Expression and metadata must contain every object cell exactly once.")
  }
  if (length(batch_col) != 1L || !batch_col %in% colnames(metadata)) {
    stop("batch_col must name an existing sample metadata column.")
  }
  batches <- as.character(metadata[cells, batch_col])
  if (anyNA(batches) || any(!nzchar(trimws(batches))) ||
      length(unique(batches)) < 2L) {
    stop("The sample column must contain at least two samples and no missing labels.")
  }
  if (is.null(features)) features <- rownames(expression)
  if (anyDuplicated(features) || !all(features %in% rownames(expression))) {
    stop("features must be unique and present in the selected expression layer.")
  }
  excluded <- features %in% exclude_features
  if (exclude_morphology) {
    excluded <- excluded |
      grepl("^(Cell|Nucleus)[_-](Area|Circularity|Solidity)$", features)
  }
  features <- features[!excluded]
  if (length(features) < 2L) stop("Fewer than two features remain.")

  # Index explicitly by cell name before transposing; AnnData uses cells as rows.
  x <- t(as.matrix(expression[features, cells, drop = FALSE]))
  storage.mode(x) <- "double"
  if (!all(is.finite(x))) stop("Selected expression contains NA, NaN, or Inf.")
  feature_sd <- apply(x, 2L, stats::sd)
  constant_features <- colnames(x)[feature_sd == 0]
  x <- x[, feature_sd > 0, drop = FALSE]
  if (ncol(x) < 2L) stop("Fewer than two variable features remain.")
  if (length(constant_features)) {
    message("Dropping constant features: ", paste(constant_features, collapse = ", "))
  }
  # Pool across samples; do not erase sample differences by centering each donor.
  if (scale_input) x <- scale(x)
  if (!all(is.finite(x))) stop("Scaling produced nonfinite measurements.")
  # R doubles can be finite yet overflow when the training matrix becomes float32.
  ranges <- apply(x, 2L, range)
  if (any(abs(ranges) > 3.4028234663852886e38)) {
    stop("Selected expression exceeds the finite float32 range; inspect scaling.")
  }
  feature_ranges <- data.frame(
    feature = colnames(x), minimum = ranges[1L, ], maximum = ranges[2L, ],
    row.names = NULL
  )
  obs <- data.frame(cytovi_batch = batches, row.names = cells)
  list(x = x, obs = obs, constant_features = constant_features,
       feature_ranges = feature_ranges)
}

#' Train CytoVI from a Seurat v5 RDS and save a new RDS with its latent embedding.
#'
#' @param input_rds Path to an existing Seurat RDS.
#' @param output_rds Path for the new RDS; must differ from input_rds.
#' @param batch_col Metadata column containing sample IDs, not biological groups.
#' @param assay Measurement assay name, default CODEX.
#' @param layer Either scale.data (used as supplied) or data (pooled z-scored).
#' @param features Optional feature names; NULL uses all features in the layer.
#' @param exclude_features Feature names to omit; defaults to character() (none).
#' @param python Optional Python executable; select it before Python initializes.
#' @param n_latent Number of CytoVI latent dimensions, default 20.
#' @param max_epochs Maximum training epochs, default 400 (the n_epochs limit).
#' @param batch_size Training minibatch size, default 256.
#' @param seed Random seed for scvi-tools, default 42. Hardware may affect results.
#' @param accelerator Lightning accelerator, e.g. auto, cpu, gpu, or mps.
#' @param model_dir Optional new directory for a reusable model plus AnnData.
#' @param overwrite Whether an existing output RDS may be replaced, default FALSE.
#' @param n_layers Hidden layers in the encoder and decoder, default 1.
#' @param dropout_rate Neural-network dropout probability in [0, 1), default 0.1.
#' @param early_stopping Whether validation performance can stop training before
#'   max_epochs, default TRUE. Does not change the separate 50-epoch KL warmup.
#' @param exclude_morphology Whether to omit Cell/Nucleus Area, Circularity, and
#'   Solidity from model input; defaults to FALSE. Excluded features remain in
#'   the saved Seurat object.
#' @param lr Optimizer learning rate, default 0.001 (CytoVI's default).
#' @param matmul_precision PyTorch float32 matrix multiplication precision:
#'   highest (default, full float32), high, or medium. The latter two can use
#'   NVIDIA Tensor Cores; these are unrelated to Google TPU accelerators.
#' @param gradient_clip_val Maximum gradient norm; default 0 disables clipping.
#' @param detect_anomaly Enable slower autograd anomaly diagnostics, default FALSE.
#' @return Invisibly, the updated Seurat object. Writes output_rds and, optionally,
#'   model_dir. Initializes Python and trains a model in the current R process.
#'   On training failure, writes output_rds.training-error.txt with runtime
#'   information, feature ranges, and the available Python error/traceback.
#' @details Requires a single exact layer containing every object cell; split
#'   layers must be joined beforehand. Does not reconstruct marker expression or
#'   change active identities. Stops before saving on invalid input, an existing
#'   CytoVI subobject, invalid latent coordinates, or dependency/training errors.
run_cytovi_rds <- function(input_rds, output_rds, batch_col,
                           assay = "CODEX", layer = "scale.data",
                           features = NULL, exclude_features = character(),
                           python = NULL, n_latent = 20L, max_epochs = 400L,
                           batch_size = 256L, seed = 42L,
                           accelerator = "auto", model_dir = NULL,
                           overwrite = FALSE, n_layers = 1L,
                           dropout_rate = 0.1, early_stopping = TRUE,
                           exclude_morphology = FALSE, lr = 0.001,
                           matmul_precision = "highest", gradient_clip_val = 0,
                           detect_anomaly = FALSE) {
  for (value in list(n_latent, n_layers, max_epochs, batch_size)) {
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
        value < 1 || value != floor(value)) {
      stop("n_latent, n_layers, max_epochs, and batch_size must be positive integers.")
    }
  }
  if (length(dropout_rate) != 1L || !is.numeric(dropout_rate) ||
      !is.finite(dropout_rate) || dropout_rate < 0 || dropout_rate >= 1) {
    stop("dropout_rate must be a number in [0, 1).")
  }
  if (length(early_stopping) != 1L || !is.logical(early_stopping) ||
      is.na(early_stopping)) {
    stop("early_stopping must be TRUE or FALSE.")
  }
  if (length(exclude_morphology) != 1L || !is.logical(exclude_morphology) ||
      is.na(exclude_morphology)) {
    stop("exclude_morphology must be TRUE or FALSE.")
  }
  if (length(lr) != 1L || !is.numeric(lr) || !is.finite(lr) || lr <= 0) {
    stop("lr must be a finite positive number.")
  }
  matmul_precision <- match.arg(matmul_precision, c("highest", "high", "medium"))
  if (length(gradient_clip_val) != 1L || !is.numeric(gradient_clip_val) ||
      !is.finite(gradient_clip_val) || gradient_clip_val < 0) {
    stop("gradient_clip_val must be a finite nonnegative number.")
  }
  if (length(detect_anomaly) != 1L || !is.logical(detect_anomaly) ||
      is.na(detect_anomaly)) {
    stop("detect_anomaly must be TRUE or FALSE.")
  }
  for (package in c("SeuratObject", "reticulate")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Install the R package '", package, "' before running this script.")
    }
  }
  if (utils::packageVersion("SeuratObject") < "5.0.0") {
    stop("SeuratObject >= 5.0.0 is required.")
  }
  layer <- match.arg(layer, c("scale.data", "data"))
  if (!file.exists(input_rds)) stop("Input RDS does not exist: ", input_rds)
  if (!dir.exists(dirname(output_rds))) stop("Create the output directory first.")
  output_path <- file.path(normalizePath(dirname(output_rds)), basename(output_rds))
  if (normalizePath(input_rds) == output_path) {
    stop("input_rds and output_rds must be different files.")
  }
  if (file.exists(output_rds) && !overwrite) stop("Output RDS already exists.")
  if (!is.null(model_dir) && file.exists(model_dir)) {
    stop("model_dir must be a new directory.")
  }
  object <- readRDS(input_rds)
  if (!inherits(object, "Seurat")) stop("Input RDS must contain a Seurat object.")
  if (!assay %in% SeuratObject::Assays(object)) stop("Assay not found: ", assay)
  if ("CytoVI" %in% names(object)) stop("A subobject named CytoVI already exists.")
  available_layers <- SeuratObject::Layers(object[[assay]], search = NA)
  if (!layer %in% available_layers) {
    stop("Exact layer '", layer, "' not found. Available: ",
         paste(available_layers, collapse = ", "),
         ". Join split layers first, or use layer = 'data' to scale transformed data.")
  }
  expression <- SeuratObject::LayerData(object, assay = assay, layer = layer)
  input <- prepare_cytovi_input(
    expression, object[[]], colnames(object), batch_col,
    features = features, exclude_features = exclude_features,
    scale_input = layer == "data", exclude_morphology = exclude_morphology
  )
  message("Training on ", nrow(input$x), " cells, ", ncol(input$x),
          " features, and ", length(unique(input$obs$cytovi_batch)), " samples.")

  if (!is.null(python)) reticulate::use_python(python, required = TRUE)
  scvi <- reticulate::import("scvi", convert = FALSE)
  external <- reticulate::import("scvi.external", convert = FALSE)
  if (!reticulate::py_has_attr(external, "CYTOVI")) {
    stop("This Python environment lacks scvi.external.CYTOVI; update scvi-tools.")
  }
  anndata <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  torch <- reticulate::import("torch", convert = FALSE)
  sys <- reticulate::import("sys", convert = FALSE)
  scvi$settings$seed <- as.integer(seed)
  torch$set_float32_matmul_precision(matmul_precision)
  runtime <- list(
    python = reticulate::py_to_r(sys$executable),
    python_version = reticulate::py_to_r(sys$version),
    scvi_version = reticulate::py_to_r(scvi$`__version__`),
    scvi_location = reticulate::py_to_r(scvi$`__file__`),
    torch_version = reticulate::py_to_r(torch$`__version__`),
    torch_location = reticulate::py_to_r(torch$`__file__`),
    cuda_available = reticulate::py_to_r(torch$cuda$is_available())
  )
  message("Python: ", runtime$python, "\n", runtime$python_version,
          "\nscvi-tools ", runtime$scvi_version, ": ", runtime$scvi_location,
          "\nPyTorch ", runtime$torch_version, ": ", runtime$torch_location,
          "\nCUDA available: ", runtime$cuda_available,
          "; float32 matmul precision: ", matmul_precision)

  # Explicit float32 C-order input avoids a read-only R-backed NumPy buffer.
  adata <- anndata$AnnData(
    X = np$array(input$x, dtype = "float32", order = "C"),
    obs = reticulate::r_to_py(input$obs)
  )
  if (!isTRUE(reticulate::py_to_r(np$isfinite(adata$X)$all()))) {
    stop("The Python float32 matrix contains NaN/Inf after conversion.")
  }
  adata$var_names <- reticulate::r_to_py(as.list(colnames(input$x)))
  external$CYTOVI$setup_anndata(adata, batch_key = "cytovi_batch")
  model <- external$CYTOVI(
    adata, n_latent = as.integer(n_latent), n_layers = as.integer(n_layers),
    dropout_rate = as.numeric(dropout_rate), protein_likelihood = "normal"
  )
  tryCatch(
    model$train(
      max_epochs = as.integer(max_epochs), batch_size = as.integer(batch_size),
      n_epochs_kl_warmup = as.integer(min(50L, max_epochs)),
      early_stopping = early_stopping, accelerator = accelerator, devices = 1L,
      lr = as.numeric(lr), precision = "32-true",
      gradient_clip_val = as.numeric(gradient_clip_val),
      detect_anomaly = detect_anomaly
    ),
    error = function(error) {
      # Capture Python's last error before calling any more Python functions.
      python_error <- reticulate::py_last_error()
      report <- c(
        conditionMessage(error), "", "Runtime:", capture.output(print(runtime)),
        "", "Input feature ranges:", capture.output(print(input$feature_ranges)),
        "", "Python error:", capture.output(print(python_error))
      )
      report_path <- paste0(output_rds, ".training-error.txt")
      tryCatch(
        writeLines(report, report_path),
        error = function(write_error) message("Could not save diagnostics: ",
                                              conditionMessage(write_error))
      )
      message("Training failed; diagnostic report path: ", report_path)
      stop(error)
    }
  )

  embedding <- as.matrix(reticulate::py_to_r(model$get_latent_representation()))
  returned_cells <- as.character(unlist(
    reticulate::py_to_r(adata$obs_names$to_list()), use.names = FALSE
  ))
  if (!identical(dim(embedding), c(length(returned_cells), as.integer(n_latent))) ||
      anyDuplicated(returned_cells) ||
      !setequal(returned_cells, colnames(object)) || !all(is.finite(embedding))) {
    stop("CytoVI returned invalid coordinates or mismatched cell names.")
  }
  rownames(embedding) <- returned_cells
  embedding <- embedding[colnames(object), , drop = FALSE]
  colnames(embedding) <- paste0("CytoVI_", seq_len(ncol(embedding)))
  object[["CytoVI"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding, assay = assay, key = "CytoVI_",
    misc = list(
      batch_col = batch_col, input_layer = layer,
      pooled_scaling = layer == "data", features = colnames(input$x),
      exclude_features = exclude_features, exclude_morphology = exclude_morphology,
      constant_features = input$constant_features,
      n_latent = n_latent, max_epochs = max_epochs, batch_size = batch_size,
      n_layers = n_layers, dropout_rate = dropout_rate,
      early_stopping = early_stopping,
      lr = lr, matmul_precision = matmul_precision, precision = "32-true",
      gradient_clip_val = gradient_clip_val, detect_anomaly = detect_anomaly,
      runtime = runtime, feature_ranges = input$feature_ranges,
      kl_warmup_epochs = min(50L, max_epochs), seed = seed,
      scvi_version = reticulate::py_to_r(scvi$`__version__`)
    )
  )
  methods::validObject(object)
  if (!is.null(model_dir)) {
    model$save(model_dir, save_anndata = TRUE, overwrite = FALSE)
  }
  saveRDS(object, output_rds)
  message("Saved Seurat object with reduction 'CytoVI': ", output_rds)
  invisible(object)
}
