#' Identify KEGG enrichment for all annotated kegg terms for an organism using Fsher tests
#' 
#' Script based on the GO enrichment analysis in k_means_figure
#' To get GO enrichment of a set of genes of interest based on universe of all expressed genes in dataset
#' Rogan Grant 2021-04-14 - 2021-12-16 
#' 
#' @param deseq_object a DESeq2 object with library normalization pre-calculated using DESeq(), vst(), etc
#' @param goi vector genes of interest in ensembl gene ID format. For fisher test: just IDs. For KS test: vector of scores with IDs as names.
#' @param kegg_organism abbreviated organism for kegg database (defaults to mouse "mmu")
#' @param gene_conversion conversion from ensembl_gene_id to external_gene_name (generally from biomart)
#' @param expression_cutoff the miminum counts detected for a given gene to be considered expressed
#' @param cores number of cores to use (defaults to 1)
#' @return a dataframe containing the significantly enriched GO terms and enrichment scores if requested
#' @export

get_kegg_enrichment = function(deseq_object, goi, kegg_organism = "mmu", 
                               gene_conversion, expression_cutoff = 1, cores = 1)
{
  #get kegg lists
  kegg_pathways = keggList("pathway", kegg_organism)
  kegg_ids = names(kegg_pathways) %>% 
    gsub("path:", "", .)
  full_kegg_names = paste(as.character(kegg_pathways), names(kegg_pathways), sep = "; ")
  
  kegg_lists = mclapply(kegg_ids, function(path){
    kegg = keggGet(path)
    genes = kegg[[1]][["GENE"]]
    if(is.null(genes))
    {
      return(NA)
    } else
    {
      #follows format ORF, description, so remove ORF and strip gene name
      genes = genes[c(F, T)] #skip ORFs
      genes = substring(genes, first = 1, last = (regexpr("\\;", genes) - 1))
    }
    
    return(genes) }, mc.cores = cores)
  names(kegg_lists) = full_kegg_names
  
  #get all expressed genes
  universe = counts(deseq_object, normalized = T) %>% 
    as.data.frame() %>% 
    rownames_to_column("ensembl_gene_id") %>% 
    left_join(., gene_conversion) %>% 
    dplyr::filter(!is.na(external_gene_name) & external_gene_name != "") %>% 
    column_to_rownames("external_gene_name") %>% 
    dplyr::select(-ensembl_gene_id) %>% 
    as.matrix() %>% 
    .[rowSums(.) > expression_cutoff, ] %>% 
    rownames()
  
  #now perform fisher tests of enrichment for each group
  enrichments = mclapply(kegg_lists, function(kegg){
    non_hits = setdiff(gene_universe, goi)
    
    #make contingency table
    hits_in_kegg = length(intersect(goi, kegg))
    hits_not_in_kegg = length(setdiff(goi, kegg))
    non_hits_in_kegg = length(intersect(non_hits, kegg))
    non_hits_not_in_kegg = length(setdiff(non_hits, kegg))
    
    contingency = rbind(c(hits_in_kegg, non_hits_in_kegg),
                        c(hits_not_in_kegg, non_hits_not_in_kegg))
    rownames(contingency) = c("in_kegg_term", "not_in_kegg_term")
    colnames(contingency) = c("DEG", "not_DEG")
    
    #get expected frequencies using chisq
    expected = suppressWarnings(chisq.test(contingency)$expected)
    
    fisher = fisher.test(contingency, alternative = "greater")
    
    #format as in topGO
    out_df = data.frame(expected = expected["in_kegg_term", "DEG"], 
                        observed = contingency["in_kegg_term", "DEG"], 
                        fold_enrichment = contingency["in_kegg_term", "DEG"] / expected["in_kegg_term", "DEG"],
                        pval = fisher$p.value) 
    return(out_df) }, mc.cores = cores)
  
  enrich_df = bind_rows(enrichments) %>% 
    dplyr::mutate(pathway = names(kegg_lists),
                  padj = p.adjust(pval, method = "fdr")) %>% 
    dplyr::arrange(padj)
  
  return(enrich_df)
}