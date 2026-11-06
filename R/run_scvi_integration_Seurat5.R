# Integrate RNA counts with scVI while retaining a complete Seurat v5 object.
# Source this file and call run_SCVI_integration(object_or_rds, ...).
# Requires R packages SeuratObject >= 5, Matrix, reticulate, and optionally renv;
# Python requires scvi-tools, anndata, numpy, pandas, scipy, and PyTorch.
# Inputs: raw RNA counts, selected variable genes, and cell-level batch labels.
# Outputs: Seurat reduction "SCVI", the returned object, and an optional new RDS.
# Split count layers are joined only in a temporary assay; source assays persist.
# This is a count-based RNA workflow, not the continuous-intensity CytoVI workflow.

#' Prepare sparse raw counts and batch labels in a common, explicit cell order.
#'
#' @param counts Gene-by-cell raw count matrix containing the selected genes.
#' @param metadata Cell metadata data.frame with cell names as row names.
#' @param cells Desired cell order, normally colnames of the complete object.
#' @param batch_col Metadata column containing batch/sample labels.
#' @return List containing sparse cells-by-genes x and aligned obs data.frame.
#' @details Stops for missing/duplicated names, invalid batch labels, nonfinite,
#'   negative or noninteger counts, or cells with zero selected-gene counts.
#'   Does not normalize, round, drop cells, or modify inputs. Materializes only
#'   the selected genes as a sparse Matrix object for transfer through reticulate.
prepare_scvi_input <- function(counts, metadata, cells, batch_col) {
  if (length(cells) < 3L || anyNA(cells) || anyDuplicated(cells) ||
      is.null(rownames(counts)) || nrow(counts) < 2L ||
      anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)) ||
      !all(cells %in% colnames(counts)) ||
      !all(cells %in% rownames(metadata))) {
    stop("Counts and metadata must contain unique names and every object cell.")
  }
  if (length(batch_col) != 1L || !batch_col %in% colnames(metadata)) {
    stop("batch_col must name an existing metadata column.")
  }
  batches <- as.character(metadata[cells, batch_col])
  if (anyNA(batches) || any(!nzchar(trimws(batches)))) {
    stop("Batch labels cannot be missing or empty.")
  }
  counts <- methods::as(counts[, cells, drop = FALSE], "dgCMatrix")
  values <- counts@x
  if (any(!is.finite(values)) || any(values < 0) ||
      any(abs(values - round(values)) > 1e-8)) {
    stop("SCVI requires finite, nonnegative integer raw RNA counts.")
  }
  if (any(Matrix::colSums(counts) == 0)) {
    stop("Some cells have zero counts across selected genes; revise features or QC.")
  }
  list(
    x = methods::as(Matrix::t(counts), "dgCMatrix"),
    obs = data.frame(scvi_batch = batches, row.names = cells)
  )
}

#' Validate and reorder scVI latent coordinates by their returned cell names.
#'
#' @param embedding Numeric cells-by-latent-dimensions matrix from scVI.
#' @param returned_cells Cell names from the corresponding AnnData object.
#' @param cells Expected Seurat cell names, in the desired output order.
#' @param n_latent Expected number of latent dimensions.
#' @return Matrix in Seurat cell order with SCVI_ dimension names.
#' @details Stops on duplicate/mismatched cells, wrong dimensions, or nonfinite
#'   coordinates. Does not alter the input matrix in the caller's environment.
align_scvi_embedding <- function(embedding, returned_cells, cells, n_latent) {
  embedding <- as.matrix(embedding)
  if (!identical(dim(embedding), c(length(returned_cells), as.integer(n_latent))) ||
      anyNA(returned_cells) || anyDuplicated(returned_cells) ||
      !setequal(returned_cells, cells) || !all(is.finite(embedding))) {
    stop("SCVI returned invalid coordinates or mismatched cell names.")
  }
  rownames(embedding) <- returned_cells
  embedding <- embedding[cells, , drop = FALSE]
  colnames(embedding) <- paste0("SCVI_", seq_len(ncol(embedding)))
  embedding
}

#' Perform scVI integration using public Seurat v5 count-layer accessors.
#'
#' @param object Seurat object or path to an RDS containing one. Counts must be
#'   raw RNA counts; variable features must already be selected unless features
#'   is provided. All object cells must be represented in the chosen assay.
#' @param python_path Python executable for reticulate; NA/NULL uses discovery.
#' @param project_path Optional renv project to load first; NA/NULL skips renv.
#' @param random_seed Random seed for scvi-tools, default 12345.
#' @param batch_col Batch/sample metadata column, default orig.ident.
#' @param n_epochs Maximum training epochs; NA/NULL uses scvi-tools' heuristic.
#' @param RDS_path Optional output RDS path; NA/NULL disables saving.
#' @param use_GPU TRUE requests one GPU, FALSE requests CPU execution.
#' @param hvgs Maximum number of existing variable genes to use, default 3000.
#'   Ignored when an explicit features vector is supplied.
#' @param n_layers Encoder/decoder hidden layers, default 1.
#' @param dropout_rate Neural-network dropout probability in [0, 1), default 0.1.
#' @param early_stopping Whether validation can stop training early, default FALSE.
#' @param n_latent Number of latent dimensions, default 10; passed to SCVI.
#' @param assay Count assay; NULL uses DefaultAssay(object). Prefer RNA when the
#'   default assay contains corrected, integrated, or non-RNA measurements.
#' @param features Optional explicit gene names, overriding variable-feature/HVG
#'   selection. An absent gene is an error rather than silently omitted.
#' @param batch_size Training and inference minibatch size, default 128.
#' @param accelerator Optional Lightning accelerator override, e.g. gpu/cpu/mps;
#'   NULL maps use_GPU to gpu or cpu. One device is used in either case.
#' @param overwrite Whether an existing SCVI subobject or output RDS can be
#'   replaced, default FALSE. An input RDS is never overwritten in place.
#' @return Invisibly, the complete Seurat object with reduction SCVI, whether or
#'   not RDS_path was supplied. Also writes RDS_path when requested.
#' @details Loads optional renv state, initializes Python, and trains scVI.
#'   Split counts/counts.* layers are joined in a temporary assay after selecting
#'   genes; the returned object's assays, identities and metadata are retained.
#'   Existing parameter names/order match the supplied legacy script. Uses
#'   accelerator/devices instead of the removed use_gpu training argument.
#'   Stops for missing dependencies, invalid counts/settings, overlapping count
#'   layers, missing features/cells, or training errors. Installs no packages.
run_SCVI_integration <- function(object,
                                 python_path = NA, project_path = NA,
                                 random_seed = 12345, batch_col = "orig.ident",
                                 n_epochs = NA, RDS_path = NA, use_GPU = FALSE,
                                 hvgs = 3000, n_layers = 1, dropout_rate = 0.1,
                                 early_stopping = FALSE, n_latent = 10,
                                 assay = NULL, features = NULL, batch_size = 128L,
                                 accelerator = NULL, overwrite = FALSE) {
  # Normalize legacy NA-valued optional arguments to NULL for Python and R use.
  if (length(python_path) == 1L && is.na(python_path)) python_path <- NULL
  if (length(project_path) == 1L && is.na(project_path)) project_path <- NULL
  if (length(n_epochs) == 1L && is.na(n_epochs)) n_epochs <- NULL
  if (length(RDS_path) == 1L && is.na(RDS_path)) RDS_path <- NULL
  integers <- list(hvgs = hvgs, n_layers = n_layers, n_latent = n_latent,
                   batch_size = batch_size)
  if (!is.null(n_epochs)) integers$n_epochs <- n_epochs
  for (name in names(integers)) {
    value <- integers[[name]]
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
        value < 1 || value != floor(value)) {
      stop(name, " must be a positive integer.")
    }
  }
  if (length(dropout_rate) != 1L || !is.numeric(dropout_rate) ||
      !is.finite(dropout_rate) || dropout_rate < 0 || dropout_rate >= 1) {
    stop("dropout_rate must be a number in [0, 1).")
  }
  for (flag in list(use_GPU, early_stopping, overwrite)) {
    if (length(flag) != 1L || !is.logical(flag) || is.na(flag)) {
      stop("use_GPU, early_stopping, and overwrite must be TRUE or FALSE.")
    }
  }
  if (!is.null(project_path)) {
    if (!requireNamespace("renv", quietly = TRUE)) stop("Install renv first.")
    renv::load(project_path)
  }
  for (package in c("SeuratObject", "Matrix", "reticulate")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Install the R package '", package, "' before running this script.")
    }
  }
  if (utils::packageVersion("SeuratObject") < "5.0.0") {
    stop("SeuratObject >= 5.0.0 is required.")
  }
  input_path <- NULL
  if (is.character(object) && length(object) == 1L) {
    input_path <- normalizePath(object, mustWork = TRUE)
    object <- readRDS(input_path)
  }
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object or RDS path.")
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)
  if (!assay %in% SeuratObject::Assays(object)) stop("Assay not found: ", assay)
  if ("SCVI" %in% names(object) && !overwrite) stop("An SCVI subobject already exists.")
  if (!is.null(RDS_path)) {
    if (!dir.exists(dirname(RDS_path))) stop("Create the output directory first.")
    output_path <- normalizePath(RDS_path, mustWork = FALSE)
    if (!is.null(input_path) && identical(input_path, output_path)) {
      stop("RDS_path must differ from the input RDS path.")
    }
    if (file.exists(RDS_path) && !overwrite) stop("Output RDS already exists.")
  }

  assay_input <- object[[assay]]
  genes <- if (is.null(features)) {
    head(SeuratObject::VariableFeatures(assay_input), as.integer(hvgs))
  } else {
    features
  }
  if (length(genes) < 2L || anyNA(genes) || anyDuplicated(genes) ||
      !all(genes %in% rownames(assay_input))) {
    stop("Select at least two existing variable genes, or supply valid features.")
  }
  # Subset only this temporary assay, leaving every assay in the result intact.
  assay_input <- subset(assay_input, features = genes)
  count_layers <- grep(
    "^counts($|[.])", SeuratObject::Layers(assay_input, search = NA), value = TRUE
  )
  if (!length(count_layers)) stop("The chosen assay has no raw counts layer.")
  count_layer <- count_layers[[1L]]
  if (length(count_layers) > 1L) {
    layer_cells <- unlist(lapply(count_layers, function(layer) {
      SeuratObject::Cells(assay_input, layer = layer)
    }), use.names = FALSE)
    if (anyDuplicated(layer_cells)) {
      stop("Raw count layers overlap in cell names; choose one count representation.")
    }
    # Use a regex to select mixed counts + counts.sample names. A distinct new
    # layer name survives JoinLayers' removal of all original source layers.
    count_layer <- "scvi_input_counts"
    while (count_layer %in% SeuratObject::Layers(assay_input, search = NA)) {
      count_layer <- paste0(count_layer, "_")
    }
    assay_input <- SeuratObject::JoinLayers(
      assay_input, layers = "^counts($|[.])", new = count_layer
    )
  }
  counts <- SeuratObject::LayerData(assay_input, layer = count_layer)
  if (!all(genes %in% rownames(counts))) stop("Selected genes are missing from counts.")
  input <- prepare_scvi_input(
    counts[genes, , drop = FALSE], object[[]], colnames(object), batch_col
  )
  message("Training scVI on ", nrow(input$x), " cells and ", ncol(input$x), " genes.")

  if (!is.null(python_path)) reticulate::use_python(python_path, required = TRUE)
  scvi <- reticulate::import("scvi", convert = FALSE)
  models <- reticulate::import("scvi.model", convert = FALSE)
  anndata <- reticulate::import("anndata", convert = FALSE)
  sparse <- reticulate::import("scipy.sparse", convert = FALSE)
  scvi$settings$seed <- as.integer(random_seed)
  # Keep the count matrix sparse across the R/Python boundary; never densify it.
  adata <- anndata$AnnData(
    X = sparse$csr_matrix(reticulate::r_to_py(input$x), dtype = "float32"),
    obs = reticulate::r_to_py(input$obs)
  )
  adata$var_names <- reticulate::r_to_py(as.list(colnames(input$x)))
  models$SCVI$setup_anndata(adata, batch_key = "scvi_batch")
  model <- models$SCVI(
    adata, n_layers = as.integer(n_layers), dropout_rate = as.numeric(dropout_rate),
    n_latent = as.integer(n_latent)
  )
  if (is.null(accelerator)) accelerator <- if (use_GPU) "gpu" else "cpu"
  model$train(
    max_epochs = if (is.null(n_epochs)) NULL else as.integer(n_epochs),
    accelerator = accelerator, devices = 1L, batch_size = as.integer(batch_size),
    early_stopping = early_stopping
  )
  returned_cells <- as.character(unlist(
    reticulate::py_to_r(adata$obs_names$to_list()), use.names = FALSE
  ))
  embedding <- align_scvi_embedding(
    reticulate::py_to_r(model$get_latent_representation(batch_size = as.integer(batch_size))),
    returned_cells, colnames(object), n_latent
  )
  object[["SCVI"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding, assay = assay, key = "SCVI_",
    misc = list(
      features = genes, count_layers = count_layers, batch_col = batch_col,
      random_seed = random_seed, n_layers = n_layers, n_latent = n_latent,
      dropout_rate = dropout_rate, n_epochs = n_epochs,
      early_stopping = early_stopping, batch_size = batch_size,
      accelerator = accelerator,
      scvi_version = reticulate::py_to_r(scvi$`__version__`)
    )
  )
  methods::validObject(object)
  if (!is.null(RDS_path)) saveRDS(object, RDS_path)
  invisible(object)
}
