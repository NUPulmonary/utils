#' Function to predict the minimum number of reads required to achieve a pre-specified read depth per cell
#' 
#' @param n_cells number of cells in the dataset
#' @param reads vector of number of reads to try
#' @param simulations number of simulations to perform
#' @param cores number of cores to use for parallel processing of simulations (defaults to 1)
#' @param random_seed self-explanatory
#' @import tidyverse
#' @return a list with "summary" containing summary stats of all simulations and "simulations" containing all individual simulation data
#' @export
determine_min_reads = function(n_cells, 
                               reads = c(miniseq_lo = 8e6, miniseq_hi = 25e6, 
                                         nextseq_P1 = 100e6, nextseq_P2 = 400e6, nextseq_P3 = 1.2e9,
                                         novaseq_SP = 800e6, novaseq_S1 = 1.6e9, novaseq_S2 = 4.1e9, novaseq_S4 = 10e9),
                               simulations = 1000, cores = 1, random_seed = 12345)
{
  library(tidyverse)
  library(parallel)
  
  #iterate through each specified cartridge and randomly assign reads to cells (labeled 1-N)
  results = lapply(reads, function(kit_reads){
    read_ids = 1:kit_reads #just name 1-N as well
    set.seed(random_seed)
    
    #iterate through desired number of simulations
    sims = mclapply(1:simulations, function(sim){
      cell_assignments = sample(1:n_cells, size = kit_reads, replace = T) # each cell given equal weight
      sim_df = tibble(read_id = read_ids, cell = cell_assignments)
      
      #fill in cells without an assigned read to account for empty
      empties = tibble(cell = setdiff(1:n_cells, sim_df$cell),
                       read_id = NA)
      
      #now zeros are accounted for
      sim_df = bind_rows(sim_df, empties) %>% 
        group_by(cell) %>% 
        dplyr::summarize(n_reads = sum(!is.na(read_id))) %>% 
        ungroup %>% 
        mutate(simulation_number = sim)
      return(sim_df) },
      mc.cores = cores)
    
    #summarize by individual simulation
    simulation_summary = bind_rows(sims) %>% 
      group_by(simulation_number) %>% 
      dplyr::summarize(median_reads = median(n_reads), #mean and median will always be reads/cells
                       mean_reads = mean(n_reads),
                       sd_reads = sd(n_reads),
                       min_reads = min(n_reads),
                       empty_cells = sum(n_reads == 0)) %>% 
      ungroup() %>% 
      mutate(kit = names(reads)[which(reads == kit_reads)],
             n_reads = kit_reads) })
  
  #finally, summarize by kit
  total_summary = bind_rows(results) %>% 
    group_by(kit) %>% 
    dplyr::summarize(median_reads = median(median_reads),
                     mean_reads = mean(mean_reads),
                     mean_sd = mean(sd_reads),
                     mean_empty_cells = mean(empty_cells),
                     global_min_reads = min(min_reads)) %>% 
    ungroup()
    
  #return all data
  out = list(summary = total_summary, simulations = results)
  return(out)
}