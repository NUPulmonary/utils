clean_sample_sheet_v2 = function(original, #path to original "raw" RNA info excel sheet
                              output_file = NA, #output path
                              reverse_complement_cols = NULL)  #names of columns to reverse complement
{
  library(tibble)
  library(dplyr)
  library(magrittr)
  require(spgs)
  
  import_raw = read_csv(original, col_names = F, skip_empty_rows = F) #import all together
  header_start = which(import_raw == "[Header]")
  data_start = which(import_raw == "[Data]")
  
  #import header as is
  header_empties = sum(is.na(import_raw[1:(data_start - 2), ])) # -2 to keep any trailing empty rows
  header = read_csv(original,
                    col_names = F,
                    skip = header_start, 
                    skip_empty_rows = F, 
                    n_max = data_start - header_empties)
  
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