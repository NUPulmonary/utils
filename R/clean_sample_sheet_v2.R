#' Function to programatically edit a bcl-convert sampleSheet to more effectively perform demultiplexing
#' 
#' @param original path to original sampleSheet to edit
#' @param output_file output path of the final, edited sampleSheet. Defaults to sampleSheet.csv in current dir
#' @param reverse_complement_cols names of columns to reverse complement
#' @param settings a list of vectors of settings (flag, value) for each setting to change. If null, runs with defaults. Note: settings from the original sampleSheet are removed.
#' @export


clean_sample_sheet_v2 = function(original, 
                              output_file = NA, 
                              reverse_complement_cols = NULL,
                              settings = NULL)  
{
  library(tibble)
  library(dplyr)
  library(magrittr)
  require(spgs)
  
  import_raw = read_csv(original, col_names = F, skip_empty_rows = F) #import all together
  header_start = which(import_raw == "[Header]")
  settings_start = which(import_raw == "[Settings]")
  data_start = which(import_raw == "[Data]")
  
  #import header as is. This will also include [reads], which we rarely if ever will need to edit.
  header_empties = sum(is.na(import_raw[1:(settings_start - 2), ])) # -2 to keep any trailing empty rows
  header = read_csv(original,
                    col_names = F,
                    skip = header_start, 
                    skip_empty_rows = F, 
                    n_max = settings_start - header_empties - 1)
  
  #reformat indexing data and add
  data = read_csv(original,
                  col_names = T,
                  skip = data_start,
                  skip_empty_rows = F)
  
  if(!is.null(reverse_complement_cols))
  {
    data = data %>% 
      dplyr::mutate(across(.cols = all_of(reverse_complement_cols), 
                           .fns = function(x) { spgs::reverseComplement(x, case = "upper") }))
  }
  
  #join together   
  add_cols = ncol(data) - ncol(header)
  for(i in 1:add_cols)
  {
    header = cbind(header, rep(NA, nrow(header)))
  }
  header = rbind(c("[Header]", rep(NA, ncol(header) - 1)), header)
  
  #generate settings chunk
  setting_chunk = rbind(c("[Settings]", rep(NA, ncol(header) - 1)))
  if(!is.null(settings))
  {
    for(i in 1:length(settings))
    {
      cur_setting = c(settings[[i]], rep(NA, ncol(header) - length(settings[[i]])))
      setting_chunk = rbind(setting_chunk, cur_setting)
    }
  }
  setting_chunk = rbind(setting_chunk, rep(NA, ncol(setting_chunk))) #add space before next chunk
  colnames(setting_chunk) = colnames(header)
  header = rbind(header, setting_chunk)
  
  header = rbind(header, c("[Data]", rep(NA, ncol(header) - 1)))
  header = rbind(header, colnames(data))
  colnames(header) = colnames(data) #necessary to bind but will get discarded later
  header = rbind(header, data)

  #export
  if(is.na(output_file)) #by default, just place in main directory
  {
    output_file = paste0(dirname(original), "/SampleSheet.csv")
  }
  write_delim(header, 
              file = output_file,
              quote = "needed",
              delim = ",",
              col_names = F,
              na = "")
}