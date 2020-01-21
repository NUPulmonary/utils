#' Run diffusion map
#'
#' Pulled from Seurat GitHub to revive for personal use
#'
#' NOTE: Prior to v2.3.4, this function used the R package diffusionMap to compute
#' the diffusion map components. This package was being archived and thus
#' RunDiffusion now uses the destiny package for the diffusion computations.
#' Please be aware that this will result in different default values as the two
#' underlying package implementations are different.
#'
#' @param object Seurat object
#' @param cells.use Which cells to analyze (default, all cells)
#' @param dims Which dimensions to use as input features
#' @param genes.use If set, run the diffusion map procedure on this subset of
#' genes (instead of running on a set of reduced dimensions). Not set (NULL) by
#' default
#' @param reduction.use Which dimensional reduction (PCA or ICA) to use for the
#' diffusion map input. Default is PCA
#' @param q.use Quantile to clip diffusion map components at. This addresses an
#' issue where 1-2 cells will have extreme values that obscure all other points.
#' 0.01 by default
#' @param max.dim Max dimension to keep from diffusion calculation
#' @param scale.clip Max/min value for scaled data. Default is 3
#' @param reduction.name dimensional reduction name, specifies the position in
#' the object$dr list. dm by default
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. DM by default
#' @param ... Additional arguments to the DiffusionMap call
#'
#' @return Returns a Seurat object with a diffusion map
#'
#' @importFrom utils installed.packages
#' @importFrom stats dist quantile
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' # Run Diffusion on variable genes
#' pbmc_small <- RunDiffusion(pbmc_small,genes.use = pbmc_small@var.genes)
#' # Run Diffusion map on first 10 PCs
#' pbmc_small <- RunDiffusion(pbmc_small,genes.use = pbmc_small@var.genes)
#' # Plot results
#' DMPlot(pbmc_small)
#' }
#'
RunDiffusion <- function(
  object,
  cells.use = NULL,
  dims = 1:5,
  genes.use = NULL,
  reduction.use = 'pca',
  q.use = 0.01,
  max.dim = 2,
  scale.clip = 10,
  reduction.name = "dm",
  reduction.key = "DM_",
  ...
) {
  # Check for destiny
  if (!'destiny' %in% rownames(x = installed.packages())) {
    stop("Please install destiny - learn more at https://bioconductor.org/packages/release/bioc/html/destiny.html")
  }
   curAssay = object@active.assay
  cells.use <- SetIfNull(x = cells.use, default = Cells(object))
  if (is.null(x = genes.use)) {
    dim.code <- object@reductions[[reduction.use]]@key
    dim.codes <- paste0(dim.code, dims)
    data.use <- FetchData(object = object, vars = dim.codes)
  } else if (! is.null(x = genes.use)) {
    genes.use <- intersect(x = genes.use, y = rownames(object@assays[[curAssay]]@scale.data))
    data.use <- MinMax(
      data = t(x = object@assays[[curAssay]]@data[genes.use, cells.use]),
      min = -1 * scale.clip,
      max = scale.clip
    )
  }
  data.dist <- as.matrix(dist(data.use))
  
  require(destiny)
  dm =  DiffusionMap(data = data.dist,
                     n_eigs = max.dim, verbose = T, ...)
  data.diffusion <- data.frame(dm@eigenvectors)
  colnames(x = data.diffusion) <- paste0(reduction.key, 1:ncol(x = data.diffusion))
  rownames(x = data.diffusion) <- cells.use
  for (i in 1:max.dim) { #quantile normalize
    x <- data.diffusion[,i]
    x <- MinMax(
      data = x,
      min = quantile(x = x, probs = q.use),
      quantile(x = x, probs = 1-q.use)
    )
    data.diffusion[, i] <- x
  }
  object@reductions[[reduction.name]] = CreateDimReducObject(
     embeddings = data.matrix(data.diffusion),
     key = reduction.key,
     assay = curAssay
  )
  return(object)
}

# Set a default value if an object is null
#
# @param x An object to set if it's null
# @param default The value to provide if x is null
#
# @return default if x is null, else x
#
SetIfNull <- function(x, default) {
  if (is.null(x = x)) {
    return(default)
  } else {
    return(x)
  }
}