# Simple R script to remove cells from Baysor output 
# which are not in both segmentation_polygons_2d.json AND segmentation.csv
# Unclear why this even happens in the first place
# Arg 1: the Baysor output director with segmentation.csv, segmentation_polygons_2d.json
# Arg 2: the path to your R environment of choice (root directory of your git repo)

#get command line arguments
args = commandArgs(TRUE)
baysor_dir = args[1]
environment = args[2]

#load environment
if(!require(renv))
{
  install.packages("renv")
}
renv::load(environment)

library(jsonlite)
library(readr)

#load files
setwd(baysor_dir)
json = fromJSON("segmentation_polygons_2d.json")
segmentation = read_csv("segmentation.csv", na = "NA")

#segmention has a prefix that has to be stripped, and unassigned have no name/cell (remove)
segmentation_prefix = unique(substring(segmentation$cell, 1, regexpr("-", segmentation$cell)))
segmentation_prefix = segmentation_prefix[segmentation_prefix != ""]
segmentation_cells = unique(substring(segmentation$cell, nchar(segmentation_prefix) + 1))
segmentation_cells = segmentation_cells[segmentation_cells != ""]

#identify cells missing from segmentation and report
ghost_cells = segmentation_cells[!(segmentation_cells %in% json$geometries$cell)]
ghost_cells = paste0(segmentation_prefix, ghost_cells)
message(paste(length(ghost_cells), "cell(s) missing from segmentation.csv:", ghost_cells))

#filter segment and re-export
if(length(ghost_cells) > 0)
{
  nrow_prev = nrow(segmentation)
  segmentation = subset(segmentation, !(cell %in% ghost_cells))
  message(paste(nrow_prev - nrow(segmentation), "rows removed."))
  write_csv(segmentation, "segmentation.csv")
}