#' Rotate a Seurat5 FOV an arbitrary number of degrees
#' Note: currently discards cell boundary info
#' 
#' @param fov a Seurat5 "FOV" object
#' @param angle_degrees rotation angle in degrees (counter-clockwise)
#' @param padding space to add between object boundary and zero point in x and y
#' @return the same FOV object, rotated
#' @import Seurat spdep dplyr
#' @export 

rotate_fov = function(fov,
                      angle_degrees,
                      padding = 100) {
  if(!("Seurat" %in% .packages()))
  {
    library(Seurat)
  }
  if(!("dplyr" %in% .packages()))
  {
    library(dplyr)
  }
  if(!("spdep" %in% .packages()))
  {
    library(spdep)
  }
  
  if(class(fov) != "FOV")
  {
    stop("FOV must be a Seurat FOV object")
  }
  
  #rotate centroids
  angle_radians = angle_degrees * pi / 180 #convert angle to radians
  new_coords = Rotation(fov$centroids@coords, angle = angle_radians)
  colnames(new_coords) = c("x", "y")
  #place back in positive coordinate space (used for molecules as well)
  x_offset = 0
  y_offset = 0
  if(min(new_coords[, "x"]) < 0)
  {
    x_offset = abs(min(new_coords[, "x"])) + padding
    new_coords[, "x"] = new_coords[, "x"] + x_offset
  }
  if(min(new_coords[, "y"]) < 0)
  {
    y_offset = abs(min(new_coords[, "y"])) + padding
    new_coords[, "y"] = new_coords[, "y"] + y_offset
  }
  rownames(new_coords) = fov$centroids@cells
  new_coords = CreateCentroids(coords = new_coords,
                               nsides = fov$centroids@nsides,
                               radius = fov$centroids@radius,
                               theta = fov$centroids@theta)
  
  #rotate molecules
  #create single dataframe and let seurat handle object creation
  new_molecules = lapply(names(fov$molecules), function(molecule){
    cur_molecule_obj = fov$molecules[[molecule]]
    new_coords = Rotation(cur_molecule_obj@coords, angle = angle_radians)
    colnames(new_coords) = c("x", "y")
    new_coords[, "x"] = new_coords[, "x"] + x_offset
    new_coords[, "y"] = new_coords[, "y"] + y_offset
    #add gene name for Seurat
    new_coords = as.data.frame(new_coords)
    new_coords$gene = molecule
    return(new_coords) })
  new_molecules = bind_rows(new_molecules)
  #cast to molecules object
  new_molecules = CreateMolecules(new_molecules)
  
  #assemble updated FOV for export
  out_fov = CreateFOV(coords = new_coords, 
                      molecules = new_molecules)
  return(out_fov)
}
