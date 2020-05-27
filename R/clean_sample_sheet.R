clean_sample_sheet = function(original, #path to original "raw" RNA info excel sheet
                              output_file = NA, #output path
                              sheet = "Basespace", #sheet name with template
                              investigator,
                              application = "NextSeq FASTQ Only", # what workflow to run
                              assay = "TruSeq HT",
                              date, #date of sequencing
                              read_length = 76, #length of sequence reads
                              species, #sample species of origin
                              adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", #sequencing adapters
                              adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
                              plate = NA, #name of "sample plate"
                              method = "RNA-seq",
                              reverse_complement_index2 = T)  #whether or not to take RC of i5 indexes
{
  require(readxl)
  require(tibble)
  require(spgs)
  import_raw = read_excel(original,
                          sheet = sheet,
                          col_names = F) #import all together
  data_start_raw = which(import_raw == "[Data]") #for reloading later
  not_empty_rows = apply(import_raw, 1, function(x){ !all(is.na(x)) })
  import_raw = import_raw[not_empty_rows, ]
  header_start = which(import_raw == "[Header]")
  data_start = which(import_raw == "[Data]")
  
  
  #reformat header
  header = import_raw[header_start:(data_start - 1), ]
  header[which(header == "ContainerID"), 1] = "Experiment Name"
  header[which(header == "FileVersion"), 1] = "IEMFileVersion"
  header[which(header == "Notes"), 1] = "Description"
  header = rbind(header, c("Investigator Name", investigator, rep(NA, ncol(header) - 2)))
  header = rbind(header, c("Workflow", "GenerateFASTQ", rep(NA, ncol(header) - 2)))
  header = rbind(header, c("Application", "NextSeq FASTQ Only", rep(NA, ncol(header) - 2)))
  header = rbind(header, c("Assay", assay, rep(NA, ncol(header) - 2)))
  header = rbind(header, c("Chemistry", "Amplicon", rep(NA, ncol(header) - 2)))
  
  #just directly add reads and settings
  header = rbind(header, rep(NA, ncol(header)))
  header = rbind(header, c("[Reads]", rep(NA, ncol(header) - 1)))
  header = rbind(header, c(read_length, rep(NA, ncol(header) - 1)))
  header = rbind(header, c(read_length, rep(NA, ncol(header) - 1)))
  header = rbind(header, rep(NA, ncol(header)))
  header = rbind(header, c("[Settings]", rep(NA, ncol(header) - 1)))
  header = rbind(header, c("Adapter", adapter1, rep(NA, ncol(header) - 2)))
  header = rbind(header, c("AdapterRead2", adapter2, rep(NA, ncol(header) - 2)))
  header = rbind(header, rep(NA, ncol(header)))
  
  #finally, reformat index data and add
  data = read_excel(original,
                    sheet = sheet,
                    col_names = T,
                    skip = (data_start_raw))
  colnames(data)[colnames(data) == "Name"] = "Sample_Name"
  colnames(data)[colnames(data) == "Index1Name"] = "I7_Index_ID"
  colnames(data)[colnames(data) == "Index1Sequence"] = "index"
  colnames(data)[colnames(data) == "Index2Name"] = "I5_Index_ID"
  colnames(data)[colnames(data) == "Index2Sequence"] = "index2"
  colnames(data)[colnames(data) == "Project"] = "Sample_Project"
  data$Sample_Plate = rep(plate, nrow(data))
  data$Description = method
  data$Species = species
  
  if(reverse_complement_index2)
  {
    data$index2 = reverseComplement(data$index2, case = "upper")
  }
  
  #join together   
  add_cols = ncol(data) - ncol(header)
  for(i in 1:add_cols)
  {
    header = cbind(header, rep(NA, nrow(header)))
  }
  header = rbind(header, c("[Data]", rep(NA, ncol(header) - 1)))
  header = rbind(header, colnames(data))
  colnames(header) = colnames(data) #necessary to bind but will get discarded later
  header = rbind(header, data)

  #export
  if(is.na(output_file)) #by default, just place in main directory
  {
    output_file = paste0(dirname(original), "/SampleSheet.csv")
  }
  write.table(x = header, 
              file = output_file,
              quote = F,
              sep = ",",
              na = "",
              row.names = F,
              col.names = F)
}