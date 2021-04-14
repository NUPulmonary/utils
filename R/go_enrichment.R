# script based on the GO enrichment analysis in k_means_figure
# To get GO enrichment of a set of genes of interest based on universe of all expressed genes in dataset   
# Rogan Grant 2021-04-14   

# Parameters
# deseq_object: a DESeq2 object with library normalization pre-calculated using DESeq(), vst(), etc
# goi: vector genes of interest in ensembl gene ID format
# go_annotations: R package with ensembl:GO annotations (defaults to mouse)

# Returns a data frame with go terms and adjusted p-values
go_enrichment = function(deseq_object,
                         goi,
                         go_annotations = "org.Mm.eg.db")
{
  require(topGO) 
  require(tidyverse)
  
  #define universe as all detected genes in dataset
  all_counts = counts(deseq_object, normalized = T)
  universe = rownames(all_counts[rowSums(all_counts) > 0, ])
  fisherTest = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")

  # 1 = selected, 0 = not selected in topGO
  selection = as.numeric(universe %in% goi)
  names(selection) = universe
  go_data = new("topGOdata", 
                ontology = "BP", 
                allGenes = selection,
                geneSel = function(x){
                  return(x == 1)},
                annot = annFUN.org, 
                mapping = go_annotations, 
                ID = "ensembl")
  
  #run Fisher test
  test_results = getSigGroups(go_data, fisherTest)
  score = as.data.frame(score(test_results))
  colnames(score) = "pval"
  score = rownames_to_column(score, var = "go_id") %>% 
    dplyr::mutate(padj = p.adjust(pval, method = "fdr")) %>% 
    dplyr::filter(padj < 0.05) %>% 
    dplyr::arrange(padj)
  
  # add descriptions using go DB
  for(i in 1:nrow(score))
  {
    score$description[i] = GOTERM[[score$go_id[i]]]@Term
  }
  
  score = score %>% 
    dplyr::mutate(full_go = paste(go_id, description))
  
  #in case of no significant go terms, return NULL
  if(nrow(score) == 0)
  {
    message("No significant enrichment detected")
    return(NULL)
  } else
  {
    return(score)
  }
}