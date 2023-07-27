#' Predict concentrations from raw MFI values using a 5PL fit
#' 
#' @param MFIs dataframe containing raw MFI values in column "MFI". All other columns are ignored.
#' @param standards dataframe containing the standard curve data with the columns "analyte", "MFI", and "expected"
#' @param return_plots toggle to output plots with standard curves and fitted points for each analyte. False by default.
#' @param plot_output_dir directory to output plots and curve summaries to. 
#' @param batch_name optional parameter to specify batch names for plot outputs
#' @param CV_cutoff cutoff for SE / predicted concentration to remove from dataset. Defaults to infinite, which turns filtering.
#' @param recode_lower_asymptote whether to assign values within the 95% CI of the lower asymptote to the lower asymptote value.
#' @param min_mfi_conc_pair a vector in the form c(desired_analyte_concentration, minimum_mfi_for_concentration) specifying the lower MFI cutoff for a given concentration at which a standard curve is considered valid.
#' @param output_summary whether to output the MFI values corresponsing to each [analyte] power of 10 for estimating cutoffs
#' @return an updated MFI dataframe with predicted values and standard errors added in "concentration" and "SE", respectively
#' @export
predict_5PL = function(MFIs, standards, return_plots = FALSE, output_dir, batch_name = NA,
                       CV_cutoff = Inf, recode_lower_asymptote = FALSE,  min_mfi_conc_pair = NULL,
                       output_summary = FALSE)
{
  library(drc)
  library(Cairo)
  library(latex2exp)
  library(sandwich)
  
  if(!dir.exists(output_dir))
  {
    dir.create(output_dir)
  }
  
  #fit curves and predict concentrations from MFI
  curves = lapply(unique(standards$analyte), function(x){
    sub = standards %>% 
      dplyr::filter(analyte == x)
    model = tryCatch({ 
      fit = drm(MFI ~ expected, data = sub,
                                 fct = LL.5(names = c("slope", "lower", "upper", "midpoint", "asymmetry")))
      #keep track of the upper limit of the 95% CI on the lower asymptote
      fit_sum = summary(fit)$coefficients
      fit$lower_asymptote = fit_sum["lower:(Intercept)", "Estimate"]
      fit$lower_asymptote_ceiling = fit_sum["lower:(Intercept)", "Estimate"] + 1.96 * fit_sum["lower:(Intercept)", "Std. Error"]
      return(fit)
    },
                     error = function(cond){
                       message(cond)
                       return(NULL) })
    return(model) })
  names(curves) = unique(standards$analyte)
  
  #remove curves that either did not converge or lack sufficient data (very rare, so far always the latter)
  bad_curves = names(curves[vapply(curves, is.null, TRUE)])
  #if specified, apply minimum MFI cutoffs for given MFI/concentration pair
  if(!is.null( min_mfi_conc_pair))
  {
    if(length( min_mfi_conc_pair) != 2)
    {
      stop(" min_mfi_conc_pair must be given in the form c(desired_analyte_concentration, minimum_mfi_for_concentration)")
    }
    below_mfi_conc_cutoff = vapply(curves[!(names(curves) %in% bad_curves)], function(curve){
      predicted_mfi = predict(curve, newdata = data.frame(conc =  min_mfi_conc_pair[1]))
      return(predicted_mfi <  min_mfi_conc_pair[2])
      }, TRUE)
    
    bad_curves = c(bad_curves, names(curves[below_mfi_conc_cutoff]))
  }
  curves = curves[!vapply(curves, is.null, TRUE)]
  
  #output powers of ten * MFI for cutoffs (if requested)
  if(output_summary == TRUE)
  {
    powers = c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6)
    predictions = lapply(powers, function(power){
      pot_df = data.frame(batch = rep(batch_name, length(curves)),
                          concentration = rep(power, each = length(curves)),
                          analyte = names(curves),
                          MFI = as.numeric(vapply(curves, function(curve){
                            predict(curve, newdata = data.frame(power = power)) }, 1))) }) %>% 
      bind_rows()
    #slightly hacky way to concatenate handle re-runs without duplication
    if(file.exists(paste0(output_dir, "/predictions_for_cutoffs.csv")))
    {
      tmp = read.csv(paste0(output_dir, "/predictions_for_cutoffs.csv"))
      predictions = bind_rows(predictions, tmp) %>% 
        unique()
    }
    write.csv(predictions, paste0(output_dir, "/predictions_for_cutoffs.csv"), row.names = F)
  }
  #also output names of bad analyte/batch combos
  discards = data.frame(curve = paste(batch_name, bad_curves))
  if(file.exists(paste0(output_dir, "/bad_curves.csv")))
  {
    tmp = read.csv(paste0(output_dir, "/bad_curves.csv"))
    discards = bind_rows(discards, tmp) %>% 
      unique()
  }
  write.csv(discards, paste0(output_dir, "/bad_curves.csv"), row.names = F)
  
  #now predict from the model
  out = MFIs %>% 
    dplyr::filter(!(analyte %in% bad_curves)) %>% 
    group_by(analyte) %>% 
    dplyr::mutate(concentration = ED(curves[[unique(analyte)]], MFI, type = "absolute", display = F,
                                     vcov. = sandwich)[,1],
                  SE = ED(curves[[unique(analyte)]], MFI, type = "absolute", display = F,
                          vcov. = sandwich)[,2],
                  CV = SE / concentration,
                  lower_asymptote_ceiling = curves[[unique(analyte)]]$lower_asymptote_ceiling,
                  lower_asymptote = curves[[unique(analyte)]]$lower_asymptote) %>% 
    #recode OOR readings
    dplyr::mutate(OOR = case_when(!is.nan(concentration) ~ "in range",
                                  # OOR high --> highest possible value
                                  is.nan(concentration) &
                                    MFI >= curves[[unique(analyte)]]$coefficients["upper:(Intercept)"] ~
                                    "high",
                                  is.nan(concentration) &
                                    MFI <= curves[[unique(analyte)]]$coefficients["lower:(Intercept)"] ~
                                    "low",
                                  TRUE ~ NA_character_),
                  concentration = case_when(recode_lower_asymptote == T & 
                                              MFI <= lower_asymptote_ceiling ~ NA_real_, #assign low-value to lower asymptote
                                            !is.nan(concentration) ~ concentration,
                                            # OOR high --> highest possible value
                                            is.nan(concentration) &
                                              MFI >= curves[[unique(analyte)]]$coefficients["upper:(Intercept)"] ~
                                              curves[[unique(analyte)]]$coefficients["upper:(Intercept)"],
                                            is.nan(concentration) &
                                              MFI <= curves[[unique(analyte)]]$coefficients["lower:(Intercept)"] ~ 0,
                                            TRUE ~ concentration)) %>%
    ungroup() %>% 
    #mark bad predictions
    dplyr::mutate(exclude = case_when(is.nan(CV) | concentration == 0 ~ FALSE, #these are OOR samples that get recoded
                                      CV > CV_cutoff ~ TRUE,
                                      CV < CV_cutoff ~ FALSE,
                                      TRUE ~ NA))
  
  if(return_plots == T)
  {
    for(ana in names(curves))
    {
      #make latex representation of curve formula
      coefficients = round(curves[[ana]]$parmMat[, 1], digits = 2)
      expression = TeX(paste("$MFI =", coefficients[2], 
                             "+ \\frac{", coefficients[3], 
                             "-", coefficients[2], 
                             "}{(1+10^{", coefficients[1], 
                             "(\\log(", ana, ")-\\log(", coefficients[4],
                             "))})^{", coefficients[5], "}}$"), output = "character")
      
                             
      plot = ggplot(NULL) +
        geom_function(fun = curves[[ana]]$curve[[1]], inherit.aes = F) +
        geom_point(data = curves[[ana]]$data[, 1:2], mapping = aes(x = expected, y = MFI), inherit.aes = F, 
                   color = "black", shape = 1) +
        geom_point(data = subset(out, analyte == ana), mapping = aes(x = concentration, y = MFI, color = exclude), 
                   inherit.aes = F, alpha = 0.7) +
        geom_errorbarh(data = subset(out, analyte == ana), 
                       mapping = aes(xmin = concentration - SE, xmax = concentration + SE, y = MFI, color = exclude),
                       alpha = 0.7) +
        geom_hline(yintercept = curves[[ana]]$lower_asymptote_ceiling, linetype = 2) +
        scale_color_manual(values = c("TRUE" = "firebrick",
                                      "FALSE" = "forestgreen")) +
        scale_x_continuous(trans = "pseudo_log", breaks = c(-1e3, -1e2, -1e1, 0, 1e1, 1e2, 1e3, 1e4)) + 
        theme_bw(base_family = "Arial") +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 16),
              axis.text.y = element_text(size = 16),
              axis.title.x = element_text(size = 24),
              axis.title.y = element_text(size = 24)) +
        annotate(geom = "text", x = 0, 
                 y = max(curves[[ana]]$data[, 2]) * 0.9, label = expression, parse = T, 
                 size = 4, hjust = 0) +
        labs(x = "Expected Concentration (pg/mL)",
             y = "MFI")
        
      
      #annotate bad plots
      if(ana %in% bad_curves)
      {
        plot = plot +
          annotate(geom = "text", x = 0, 
                   y = max(curves[[ana]]$data[, 2]) * 0.8, label = "EXCLUDED", 
                   size = 4, hjust = 0)
      }
      
      outname = paste0(output_dir, "/", batch_name, "_", gsub("\\.|\\/", "-", ana), ".pdf")
      CairoPDF(file = outname, width = 8, height = 8)
      plot(plot)
      dev.off()
      
    }
  }
  
  #now perform exclusions after plotting full dataset
  out = out %>% 
    dplyr::filter(exclude == FALSE)
       
  return(out)
}

#' Load wide-form data from EveTech
#' 
#' Function to read data output from EveTech, predict concentration using 5PL model, 
#' and join into a single usable dataset for further analysis in R or for CSV output
#' 
#' @param datafiles a named vector of paths to evetech data output in xslx format, where the name is the batch ID
#' @param sample_metadata dataframe or path to dataframe with minimal columns: tube location (e.g. A1), sample type (BAL or serum), and BAL ID (tc_pt_study_id)
#' @param patient_metadata the current distribution of the script metadata from MS
#' @param bind_by_sample a vector of column names to bind by for sample metadata a la *_join(). Defaults to NULL, which is all matching columns
#' @param bind_by_metadata a vector of column names to bind by for patient metadata a la *_join(). Defaults to NULL, which is all matching columns
#' @param analyte_conv a conversion table of analyte names to common names (optional)
#' @param bind_by_metadata a vector of column names to bind by for analyte conversions a la *_join(). Defaults to NULL, which is all matching columns
#' @param collapse_replicates whether or not to summarize techical replicates by mean and SD
#' @param return_plots toggle to output plots with standard curves and fitted points for each analyte. False by default.
#' @param output_dir directory to output plots and curve summaries to.
#' @param CV_cutoff cutoff for SE / predicted concentration to remove from dataset. Defaults to infinite, which turns filtering.
#' @param recode_lower_asymptote whether to assign values within the 95% CI of the lower asymptote to the lower asymptote value
#' @param  min_mfi_conc_pair a vector in the form c(desired_analyte_concentration, minimum_mfi_for_concentration) specifying the lower MFI cutoff for a given concentration at which a standard curve is considered valid.
#' @param output_summary whether to output the MFI values corresponsing to each [analyte] power of 10 for estimating cutoffs
#' @return a single dataframe of concentrations and metdata in long format
#' @export
read_evetech_predict = function(datafiles, sample_metadata = NULL, patient_metadata = NULL, 
                        bind_by_sample = NULL, bind_by_metadata = NULL,
                        analyte_conv = NULL, bind_by_analyte = NULL,
                        collapse_replicates = T, return_plots = F, output_dir,
                        CV_cutoff = Inf, recode_lower_asymptote = FALSE,  min_mfi_conc_pair = NULL,
                        output_summary = FALSE)
{
  library(tidyverse)
  library(readxl)
  
  # load files and join together
  # determine NA strings based on params
  nas = c("", "na", "NA", "N/A", "***", "---")
  data_list = vector(length = length(datafiles), mode = "list") #preallocate for speed
  for(i in 1:length(datafiles))
  {
    #determine column names on the fly due to bad formatting
    colnames = read_excel(path = datafiles[i], skip = 7, n_max = 1, sheet = "Obs Conc") %>% 
      dplyr::rename(Type = ...1, Well = ...2, Description = ...3) %>% 
      colnames() %>% 
      #remove numbering
      gsub(" \\(\\d+\\)", "", .)

    raw_mfi = suppressMessages(read_excel(path = datafiles[i], na = nas, skip = 8, col_names = colnames, sheet = "FI")) %>%
      #subset to just experimental
      dplyr::filter(grepl("X\\d+", Type) & !duplicated(.)) %>% #all of their data is duplicated for inexplicable reasons
      pivot_longer(cols = 4:ncol(.), names_to = "analyte", values_to = "MFI") %>% 
      dplyr::mutate(batch = names(datafiles[i]),
                    sample_origin = toupper(substring(names(datafiles[i]), regexpr("_", names(datafiles[i])) + 1)),
                    Type = toupper(Type),
                    Well = toupper(Well),
                    MFI = as.numeric(MFI),
                    Description = toupper(Description))
    
    standard_mfi = suppressMessages(read_excel(path = datafiles[i], na = nas, skip = 8, col_names = colnames, sheet = "FI")) %>%
      #subset to just experimental
      dplyr::filter(grepl("eS\\d+", Type) & !duplicated(.) & !grepl(",", Well)) %>% #mean columns are grouped with a comma. Ignore to avoid resampling.
      pivot_longer(cols = 4:ncol(.), names_to = "analyte", values_to = "MFI") %>% 
      dplyr::mutate(batch = names(datafiles[i]),
                    sample_origin = toupper(substring(names(datafiles[i]), regexpr("_", names(datafiles[i])) + 1)),
                    Type = toupper(Type),
                    MFI = as.numeric(MFI),
                    Well = toupper(Well),
                    Description = toupper(Description))
    standard_expected = suppressMessages(read_excel(path = datafiles[i], na = nas, skip = 8, col_names = colnames, sheet = "Exp Conc")) %>%
      #subset to just experimental
      dplyr::filter(grepl("eS\\d+", Type) & !duplicated(.) & !grepl(",", Well)) %>% 
      pivot_longer(cols = 4:ncol(.), names_to = "analyte", values_to = "expected") %>% 
      dplyr::mutate(batch = names(datafiles[i]),
                    sample_origin = toupper(substring(names(datafiles[i]), regexpr("_", names(datafiles[i])) + 1)),
                    Type = toupper(Type),
                    Well = toupper(Well),
                    expected = as.numeric(expected),
                    Description = toupper(Description))
    standards = left_join(standard_mfi, standard_expected)
    rm(standard_mfi, standard_expected)
    
    #fit curves and predict concentrations from MFI
    raw_mfi = predict_5PL(MFIs = raw_mfi, standards = standards,
                          return_plots = return_plots, output_dir = output_dir,
                          batch_name = names(datafiles[i]), CV_cutoff = CV_cutoff,
                          recode_lower_asymptote = recode_lower_asymptote,
                           min_mfi_conc_pair =   min_mfi_conc_pair, output_summary = output_summary)
    
    #add sample metadata?
    if(!is.null(sample_metadata))
    {
      raw_mfi = left_join(raw_mfi, sample_metadata, by = bind_by_sample, na_matches = "never")
    }
    
    #add clinical metadata?
    if(!is.null(patient_metadata))
    {
      raw_mfi = left_join(raw_mfi, patient_metadata, by = bind_by_metadata, na_matches = "never")
    }
    
    #add analyte common names?
    if(!is.null(analyte_conv))
    {
      display_name_col = setdiff(colnames(analyte_conv), colnames(raw_mfi))
      analyte_name_col = setdiff(colnames(analyte_conv), display_name_col)
      raw_mfi = left_join(raw_mfi, analyte_conv, by = bind_by_analyte, na_matches = "never") %>%
        #take care of missing analytes, if necessary
        dplyr::mutate(!!display_name_col := ifelse(is.na(get(display_name_col)),
                                                   yes = get(analyte_name_col),
                                                   no = get(display_name_col)))
    }
    
    data_list[[i]] = raw_mfi
  }
  
  #get rid of any empty dfs
  keepers = which(vapply(data_list, nrow, 1) > 0)
  data_list = data_list[keepers]
  
  #finally join them all together
  out = bind_rows(data_list) %>% 
    dplyr::mutate(unique_id = paste(tc_pt_study_id, sample_origin, sep = "_"),
                  sample_origin = factor(sample_origin),
                  originator = "Eve Technologies") %>% 
    dplyr::select(-matches("case_number"))
  
  #test for duplication, report, and take means
  dup_test = out %>% 
    dplyr::select(tc_pt_study_id, sample_origin, analyte)
  if(any(duplicated(dup_test)))
  {
    warning(paste("Warning:", sum(duplicated(dup_test)), "duplicate entries found"))
  }
  
  if(collapse_replicates == T)
  {
    #summary does nothing if there are no duplicates, but keep for consistency
    out = out %>% 
      group_by(across(-c(concentration, SE, OOR, MFI, Well,  batch, 
                         lower_asymptote_ceiling, lower_asymptote))) %>% 
      dplyr::summarize(mean_concentration = mean(concentration, na.rm = T),
                       sd_concentration = sd(concentration, na.rm = T),
                       cv_concentration = sd(concentration, na.rm = T) / mean(concentration, na.rm = T),
                       individual_vals = list(concentration),
                       OORs = list(OOR),
                       individual_SEs = list(SE),
                       corresponding_wells = list(Well),
                       n_replicates = n()) %>% 
      ungroup()
  }
  
  out = out %>% dplyr::select(-analyte)
  return(out)
}

#' Load CSV output from in-house Luminex assays   
#' 
#' Function to read data output from EveTech, predict concentration using 5PL model, 
#' and join into a single usable dataset for further analysis in R or for CSV output
#' 
#' @param datafiles a named vector of paths to evetech data output in xslx format, where the name is the batch ID
#' @param patient_metadata the current distribution of the script metadata from MS
#' @param bind_by_metadata a vector of column names to bind by for patient metadata a la *_join(). Defaults to NULL, which is all matching columns
#' @param analyte_conv a conversion table of analyte names to common names (optional)
#' @param bind_by_metadata a vector of column names to bind by for analyte conversions a la *_join(). Defaults to NULL, which is all matching columns
#' @param collapse_replicates whether or not to summarize techical replicates by mean and SD
#' @param min_events the minimum number of beads to keep a given analyte. Defaults to 50.   
#' @param correct_by_volume whether or not to adjust concentrations based on sample volume (divide by volume/25)
#' @param return_plots toggle to output plots with standard curves and fitted points for each analyte. False by default.
#' @param output_dir directory to output plots and curve summaries to.
#' @param CV_cutoff cutoff for SE / predicted concentration to remove from dataset. Defaults to infinite, which turns filtering.
#' @param recode_lower_asymptote whether to assign values within the 95% CI of the lower asymptote to the lower asymptote value
#' @param  min_mfi_conc_pair a vector in the form c(desired_analyte_concentration, minimum_mfi_for_concentration) specifying the lower MFI cutoff for a given concentration at which a standard curve is considered valid.
#' @param output_summary whether to output the MFI values corresponsing to each [analyte] power of 10 for estimating cutoffs
#' @return a single dataframe of concentrations and metdata in long format
#' @export
read_inhouse_predict = function(data_dirs, patient_metadata = NULL, 
                                bind_by_metadata = NULL, analyte_conv = NULL, 
                                bind_by_analyte = NULL, collapse_replicates = T, 
                                min_events = 50, correct_by_volume = T,
                                return_plots = F, output_dir,
                                CV_cutoff = Inf, recode_lower_asymptote = FALSE,
                                 min_mfi_conc_pair = NULL, output_summary = FALSE)
{
  library(tidyverse)
  library(readxl)
  
  #need this to read massive datafiles (should cause no issues)
  Sys.setenv("VROOM_CONNECTION_SIZE" = 1e6)
  
  # load files and join together
  data_list = vector(length = length(data_dirs), mode = "list") #preallocate for speed
  for(i in 1:length(data_dirs))
  {
    #get all relevant assay outputs   
    #all located in ./Plate[12]_PX\\d{2}
    datafiles = list.files(path = data_dirs[[i]], 
                           recursive = T,
                           full.names = T,
                           pattern = "^\\d{8}.+\\.csv$")
    
    #get MFI data
    raw_mfis = lapply(datafiles, function(x){
      plate_name = dirname(x) %>% 
        basename(.)
      
      #find starts and ends of relevant sections
      mfi_start = suppressWarnings(readLines(x)) %>% 
        grepl("\"DataType:\",\"Net MFI\"|DataType:,Net MFI,", .) %>% 
        which(.) + 1
      mfi_end = suppressWarnings(readLines(x)) %>% 
        grepl("\"DataType:\",\"Count\"|DataType:,Count", .) %>% 
        which(.) - 2
      
      #keep only samples, not standards for this purpose
      out = read_csv(x, 
                     skip = (mfi_start - 1), 
                     n_max = (mfi_end - mfi_start), 
                     show_col_types = F) %>% 
        dplyr::filter(grepl("Unknown", Sample)) %>% 
        dplyr::mutate(plate = plate_name) %>% 
        #perform bead number filtering
        dplyr::filter(`Total Events` >= min_events) %>% 
        dplyr::select(-`Total Events`) %>% 
        pivot_longer(cols = -c(Location, Sample, plate),
                     names_to = "analyte",
                     values_to = "MFI") 
      
      return(out)
    } )
    
    raw_mfi = bind_rows(raw_mfis) %>% 
      dplyr::mutate(Sample = factor(Sample),
                    plate = factor(plate),
                    analyte = gsub("^\\d+ ", "", analyte), #handle numbered analytes from new software
                    batch = basename(data_dirs[[i]]))
    rm(raw_mfis)
    
    #collect standard curve data (MFI and expected)
    standards = lapply(datafiles, function(x){
      plate_name = dirname(x) %>% 
        basename(.)
      
      #find starts and ends of relevant sections
      mfi_start = suppressWarnings(readLines(x)) %>% 
        grepl("\"DataType:\",\"Net MFI\"|DataType:,Net MFI,", .) %>% 
        which(.) + 1
      mfi_end = suppressWarnings(readLines(x)) %>% 
        grepl("\"DataType:\",\"Count\"|DataType:,Count", .) %>% 
        which(.) - 2
      expected_start = suppressWarnings(readLines(x)) %>% 
        grepl("\"DataType:\",\"Standard Expected Concentration\"|DataType:,Standard Expected Concentration", .) %>% 
        which(.) + 1
      expected_end = suppressWarnings(readLines(x)) %>% 
        grepl("\"DataType:\",\"Control Expected Concentration\"|DataType:,Standard Expected Concentration", .) %>% 
        which(.) - 2
      
      #keep only samples, not standards for this purpose
      MFI = read_csv(x, 
                     skip = (mfi_start - 1), 
                     n_max = (mfi_end - mfi_start),
                     show_col_types = F) %>% 
        dplyr::filter(grepl("Standard", Sample)) %>% 
        dplyr::mutate(plate = plate_name) %>% 
        #perform bead number filtering
        dplyr::filter(`Total Events` >= min_events) %>% 
        dplyr::select(-`Total Events`) %>% 
        pivot_longer(cols = -c(Location, Sample, plate),
                     names_to = "analyte",
                     values_to = "MFI") 
      
      expected = read_csv(x, 
                          skip = (expected_start - 1), 
                          n_max = (expected_end - expected_start),
                          show_col_types = F) %>% 
        dplyr::rename(Sample = Reagent) %>% 
        dplyr::filter(grepl("Standard", Sample)) %>% 
        dplyr::mutate(plate = plate_name) %>% 
        pivot_longer(cols = -c(Sample, plate),
                     names_to = "analyte",
                     values_to = "expected") %>% 
        dplyr::mutate(expected = as.double(expected))
      
      out = left_join(MFI, expected)
      
      return(out)
    } )
    
    standard = bind_rows(standards) %>% 
      dplyr::mutate(Sample = factor(Sample),
                    plate = factor(plate),
                    analyte = gsub("^\\d+ ", "", analyte), #handle numbered analytes from new software
                    analyte = factor(analyte))
    rm(standards)
    
    #fit curves and predict concentrations from MFI
    raw_mfi = predict_5PL(MFIs = raw_mfi, standards = standard,
                          return_plots = return_plots, output_dir = output_dir,
                          batch_name = basename(data_dirs[[i]]), CV_cutoff = CV_cutoff,
                          recode_lower_asymptote = recode_lower_asymptote, 
                           min_mfi_conc_pair =   min_mfi_conc_pair, output_summary = output_summary)
    
    #add sample metadata
    sample_metadata = read_csv(paste0(data_dirs[[i]], "/Batch_definition.csv"), show_col_types = F) %>% 
      dplyr::select(tc_pt_study_id, sample_origin = `Sample type`, coordinate, sample_volume = `Sample volume (ul)`) 
    raw_mfi = raw_mfi %>% 
      dplyr::mutate(coordinate = substring(Location, regexpr("[[:upper:]]", Location), nchar(Location) - 1)) %>% 
      #no definition = bad sample
      right_join(., sample_metadata)
    
    #add clinical metadata?
    if(!is.null(patient_metadata))
    {
      raw_mfi = left_join(raw_mfi, patient_metadata, by = bind_by_metadata, na_matches = "never")
    }
    
    #add analyte common names?
    if(!is.null(analyte_conv))
    {
      display_name_col = setdiff(colnames(analyte_conv), colnames(raw_mfi))
      analyte_name_col = setdiff(colnames(analyte_conv), display_name_col)
      raw_mfi = left_join(raw_mfi, analyte_conv, by = bind_by_analyte, na_matches = "never") %>%
        #take care of missing analytes, if necessary
        dplyr::mutate(!!display_name_col := ifelse(is.na(get(display_name_col)),
                                    yes = get(analyte_name_col),
                                    no = get(display_name_col)))
    }
    
    data_list[[i]] = raw_mfi
  }
  
  #get rid of any empty dfs
  keepers = which(vapply(data_list, nrow, 1) > 0)
  data_list = data_list[keepers]
  
  #finally join them all together
  out = bind_rows(data_list) %>% 
    dplyr::mutate(unique_id = paste(tc_pt_study_id, sample_origin, sep = "_"),
                  sample_origin = factor(sample_origin),
                  originator = "Northwestern University")
  
  #correct concentrations by final volume?
  if(correct_by_volume == T)
  {
    out = out %>% 
      dplyr::mutate(concentration = concentration / (sample_volume / 25))
  }
  
  #test for duplication, report, and take means
  dup_test = out %>% 
    dplyr::select(tc_pt_study_id, sample_origin, analyte)
  if(any(duplicated(dup_test)))
  {
    warning(paste("Warning:", sum(duplicated(dup_test)), "duplicate entries found"))
  }
  
  if(collapse_replicates == T)
  {
    #summary does nothing if there are no duplicates, but keep for consistency
    out = out %>% 
      dplyr::rename(Well = coordinate) %>% 
      group_by(across(-c(concentration, SE, OOR, MFI, batch, Well, Location, Sample, 
                         plate, sample_volume, lower_asymptote_ceiling, lower_asymptote))) %>% 
      dplyr::summarize(mean_concentration = mean(concentration, na.rm = T),
                       sd_concentration = sd(concentration, na.rm = T),
                       cv_concentration = sd(concentration, na.rm = T) / mean(concentration, na.rm = T),
                       individual_vals = list(concentration),
                       individual_SEs = list(SE),
                       OORs = list(OOR),
                       corresponding_wells = list(Well),
                       n_replicates = n()) %>% 
      ungroup()
  }
  
  out = out %>% dplyr::select(-analyte)
  return(out)
}