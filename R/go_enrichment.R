#' Identify GO enrichment of a list of gene hits
#' 
#' Script based on the GO enrichment analysis in k_means_figure
#' To get GO enrichment of a set of genes of interest based on universe of all expressed genes in dataset
#' Rogan Grant 2021-04-14 - 2021-12-16 
#' 
#' @param deseq_object a DESeq2 object with library normalization pre-calculated using DESeq(), vst(), etc
#' @param goi vector genes of interest in ensembl gene ID format
#' @param go_annotations R package with ensembl:GO annotations (defaults to mouse)
#' @param return_fold_enrichment whether or not to return enrichment data for each GO hit. Defaults to FALSE for backwards compatibility.
#' @param expression_cutoff the miminum counts detected for a given gene to be considered expressed
#' @return a dataframe containing the significantly enriched GO terms and enrichment scores if requested
#' @export

go_enrichment = function(deseq_object,
                         goi,
                         go_annotations = "org.Mm.eg.db",
                         return_fold_enrichment = FALSE,
                         expression_cutoff = 1)
{
  library(topGO) 
  library(tidyverse)
  
  #define universe as all detected genes in dataset
  all_counts = counts(deseq_object, normalized = T)
  universe = rownames(all_counts[rowSums(all_counts) >= expression_cutoff, ])
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
  if(return_fold_enrichment == TRUE)
  {
    score = GenTable(go_data, 
                     pval = test_results, 
                     orderBy = "pval", 
                     topNodes = length(test_results@score)) %>%  #generally just want all filtered terms; Inf returns error
      dplyr::rename(go_id = GO.ID,
                    description = Term) %>% 
      mutate(padj = p.adjust(pval, method = "fdr"),
             fold_enrichment = Significant / Expected,
             term_coverage = Significant / Annotated,
             full_go = paste(go_id, description)) %>% #
      dplyr::filter(padj < 0.05)
  } else
  {
    score = as.data.frame(score(test_results))
    colnames(score) = "pval"
    score = rownames_to_column(score, var = "go_id") %>% 
      dplyr::mutate(padj = p.adjust(pval, method = "fdr")) %>% 
      dplyr::filter(padj < 0.05) %>% 
      dplyr::arrange(padj)
  }
  
  #in case of no significant go terms, return NULL
  if(nrow(score) == 0)
  {
    message("No significant enrichment detected")
    return(NULL)
  } else
  {
    # add descriptions using go DB
    for(i in 1:nrow(score))
    {
      score$description[i] = GOTERM[[score$go_id[i]]]@Term
    }
    
    score = score %>% 
      dplyr::mutate(full_go = paste(go_id, description))
    return(score)
  }
}