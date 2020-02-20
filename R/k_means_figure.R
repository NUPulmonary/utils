#function to automatically generate an optimized k means plot
#optionally, with GO terms for each cluster
# dge: a DESeq2 dataset with DEA already run
# qval_cutoff: maximum q-value to be considered significant
# genes of interest: a vector of genes to consider (overrides other options)
   
k_means_figure = function(dge,
                          qval_cutoff = 0.05,
                          genes_of_interest,
                          