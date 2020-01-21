# generalized function to calculate RNA velocity for a Seurat dataset
# returns the same seurat object with an rna velocity slot added

# object: seurat object
# loomPaths: Vector of full paths to the velocyto output loom files from the python version of velocyto
# assay: which assay to use (e.g. SCT); current defaults to DefaultAssay
# reduction: which embeddings to use (e.g. tSNE); currently defaults to UMAP

getRNAVelocity = function(object, loomPaths, reduction = "umap", cores = 1)
{
  #dependencies
  require(Seurat)
  require(velocyto.R)
  require(future)
  require(doParallel)
  registerDoParallel(cores)
  Sys.setenv(MC_CORES=cores)
  require(dplyr)
  require(tibble)
   
  #load into velocyto
  looms = mclapply(loomPaths, read.loom.matrices) 
  names(looms) = gsub("\\.loom", "", basename(loomPaths))
  
  #make spliced matrix
  emats = mclapply(looms, function(x){
    emat = x$spliced
    colnames(emat) = substr(gsub("\\:", "___", gsub("\\_", "", colnames(emat))), 1, 23) #put in Seurat format
    return(emat) })
  emat = do.call(cbind, emats)
  emat = emat[, colnames(emat) %in% Cells(object)] #remove cells that have previously been discarded in Seurat
  rm(emats)
  
  #make unpliced matrix
  nmats = mclapply(looms, function(x){
    nmat = x$unspliced
    colnames(nmat) = substr(gsub("\\:", "___", gsub("\\_", "", colnames(nmat))), 1, 23) #put in Seurat format
    return(nmat) })
  nmat = do.call(cbind, nmats)
  nmat = nmat[, colnames(nmat) %in% Cells(object)] #remove cells that have previously been discarded in Seurat
  rm(nmats)
  rm(looms)
  gc()
  
  #recommended by Peter Kharchenko
  clusterLabels = object@meta.data$seurat_clusters
  names(clusterLabels) = rownames(object@meta.data)
  emat = filter.genes.by.cluster.expression(emat, clusterLabels, min.max.cluster.average = 0.5)
  nmat = filter.genes.by.cluster.expression(nmat, clusterLabels, min.max.cluster.average = 0.05)
  
  #get RNA velocity
  reductionDistances = as.dist(1 - cor(t(object@reductions[[reduction]]@cell.embeddings)))
  velocity = gene.relative.velocity.estimates(emat,
                                              nmat,
                                              cell.dist = reductionDistances,
                                              fit.quantile = 0.02, #this is not well explained in the vignette
                                              n.cores=cores)
  velocityName = paste0("velocity_", reduction)
  object@reductions[[velocityName]] = velocity
  
  return(object)
}


#function for plotting RNA velocity data after calculation with getRNAVelocity
plotVelocity = function(object, reduction = "umap", dims = c(1,2), group.by = "ident", cores = 1)
{
  require(ggplot2)
   embeds = object@reductions[[reduction]]@cell.embeddings
   seuratPlot = DimPlot(object = object, dims = dims, reduction = reduction, group.by = "ident")
   plotData = ggplot_build(seuratPlot)$data
   colors = as.list(plotData[[1]]$colour)
   names(colors) = rownames(embeds)
   
   velocity = object@reductions[[paste0("velocity_", reduction)]]
   out = show.velocity.on.embedding.cor(embeds, velocity,
                                        cell.colors=ac(colors,alpha=0.5),
                                        cex=0.8,arrow.scale=2,show.grid.flow=T,
                                        min.grid.cell.mass=1.0,grid.n=50,arrow.lwd=1,
                                        do.par=F,cell.border.alpha = 0.1,
                                        n.cores=cores,main="Cell Velocity")
   return(out)
}
   