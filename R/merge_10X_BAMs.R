#' SCRIPT to merge cellranger BAMs from multiple runs into single BAM files
#' Useful for merging repeated samples such as with miniseq pre-sequencing
#' 
#' @param input_dirs vector of top-level directories containing the cellranger output subdirectories
#' @param samples_of_interest vector of sample names to combine, e.g. SC001. If NULL, takes all from first directory.
#' @param output_dir directory to output merged BAM files
#' @param identity_table a dataframe with ncol = length(input_dirs) and nrow = length(samples of interest). Rownames denote final library name, rows are vectors of equivalent samples. If NULL, done automatically by name
#' @param overwrite_files whether or not to overwrite destination BAM files. Defaults to TRUE 
#' @export

merge_10X_BAMs = function(input_dirs,
                          samples_of_interest = NULL,
                          output_dir,
                          identity_table = NULL,
                          overwrite_files = TRUE)
{
  require(Rsamtools)
  
  #create output directory if necessary
  if(!dir.exists(output_dir))
  {
    dir.create(output_dir, recursive = T)
  }
  
  #if samples are not specified, take all from first directory
  if(is.null(samples_of_interest))
  {
    samples_of_interest = basename(list.dirs(input_dirs[1], recursive = F))
  }
  
  #if identity table is not specified, assume all samples match in all directories
  if(is.null(identity_table))
  {
    identity_table = matrix(nrow = length(samples_of_interest),
                            ncol = length(input_dirs))
    for(i in 1:length(samples_of_interest))
    {
      sample = samples_of_interest[i]
      cur = paste(input_dirs, sample, sep = "/")
      identity_table[i, ] = cur
    }
    rownames(identity_table) = samples_of_interest
  }
  
  #ensure all of the paths are valid
  all_files = as.vector(identity_table)
  valid_vector = vapply(all_files, dir.exists, FUN.VALUE = TRUE)
  if(any(valid_vector == FALSE))
  {
    error_message = paste0("following pilepaths are invalid\n", 
                          as.character(paste(all_files[!valid_vector], collapse = "\n")))
    stop(error_message)
  }
  
  #if all is well, merge the bam files
  for(i in 1:nrow(identity_table))
  {
    files_vector = identity_table[i, ]
    sample = rownames(identity_table)[i]
    mergeBam(files = files_vector,
             destination = paste0(output_dir, "/", sample, ".bam"))
  }
}