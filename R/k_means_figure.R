#function to automatically generate an optimized k means plot
#optionally, with GO terms for each cluster
# dge: a DESeq2 dataset with DEA already run
# qval_cutoff: maximum q-value to be considered significant
# genes_of_interest: a vector of genes to consider (overrides other options)
# display_go_terms: whether or not to display go terms for each cluster
# max_go_terms: maximum number of significant GO terms to display for each cluster
# design: a character representation of the design for an ANOVA-like test (LRT)
# cores: number of cores to run in parallel
# max_k: maximum value of k to consider
# actual k to use in final plot
# colnames: whether or not to display column labels in heatmap
# legend_factors: vector of factors to add to heatmap legend (must be in des metadata)
# tidy_go: if true, join go terms into a tidy data frame

#for extracting counts for genes of interest from a DESeqDataSet
construct_goi_matrix = function(dge,
                                qval_cutoff = 0.05,
                                genes_of_interest = NULL,
                                design = NULL,
                                cores = 1,
                                fitType = "local",
                                minReps = 7,
                                baseMeanCutoff = 0)
{
  
 if(cores > 1)
 {
  library(BiocParallel)
   BiocParallel::register(MulticoreParam(cores))
 }
  clusterAnnos = NULL #updated later if necessary
  
  counts_mat = DESeq2::counts(dge, normalized = T)
  if(is.na(qval_cutoff) || (is.null(design) && is.null(genes_of_interest))) #all cases to include all genes
  {
    qval_cutoff = 1
  }
  
  #determine subset of genes to use
  if(is.null(genes_of_interest) && qval_cutoff < 1)
  {
    #perform anova-like test to identify genes which vary significantly across factor(s) of interest
    design = as.formula(design)
    DESeq2::design(dge) = design
    deseq_results = DESeq(dge, 
                          test = "LRT", 
                          reduced = ~ 1,
                          parallel = T,
                          fitType = fitType, 
                          minReplicatesForReplace = minReps)
    deseq_results = as.data.frame(results(deseq_results, alpha = qval_cutoff, parallel = TRUE))
    genes_of_interest = rownames(subset(deseq_results, padj < qval_cutoff), baseMean > baseMeanCutoff)
  } else if(is.null(genes_of_interest))
  {
    genes_of_interest = rownames(counts_mat)
  }
  
  #subset counts down to chosen genes
  counts_mat = counts_mat[genes_of_interest, ]
  
  #z-score matrix (necessary for k means)
  counts_mat = t(base::scale(t(counts_mat))) #transpose issue is annoying
  
  return(counts_mat)
}
  

#for heuristically determining an ideal k (run this first)
k_elbow = function(dge,
                   qval_cutoff = 0.05,
                   genes_of_interest = NULL,
                   design = NA,
                   cores = 1,
                   random_seed = 12345,
                   max_k = 50,
                   minReps = 7) #this is the default for DESeq2
{
  library(ggplot2)
  library(tidyverse)
  options(gsubfn.engine = "R")
  
  counts_mat = construct_goi_matrix(dge = dge,
                                    qval_cutoff = qval_cutoff,
                                    genes_of_interest = genes_of_interest,
                                    design = design,
                                    cores = cores,
                                    minReps = minReps)
  
  #now run kmeans for all values of k, and find sums of squared differences within each cluster for each k
  sums_of_squares = mclapply(1:max_k, function(k){
    set.seed(random_seed)
    kmeans_result = kmeans(x = counts_mat, 
                           centers = k, 
                           iter.max = 1000,
                           nstart = 25)
    return(kmeans_result$tot.withinss)}, mc.cores = cores)
  sums_of_squares = unlist(sums_of_squares)
  ss_df = data.frame(k = 1:max_k, total_within_ss = sums_of_squares)
  
  #now plot and return
  plot = ggplot(ss_df, aes(x = k, y = total_within_ss)) +
    geom_point() +
    xlab("Clusters (k)") +
    ylab("Total within-cluster sum-of-squares")
  
  return(plot)
}
      
k_means_figure = function(dge,
                          qval_cutoff = 0.05,
                          genes_of_interest = NULL,
                          design = NA,
                          cores = 1,
                          k,
                          display_go_terms = T,
                          return_go_terms = F,
                          max_go_terms = 5,
                          colnames = F,
                          legend_factors = NULL,
                          go_annotations = "org.Mm.eg.db",
                          go_ontology = "BP",
                          ensembl_db = "mmusculus_gene_ensembl",
                          cluster_columns = T,
                          return_genes = F,
                          label_fontsize = 6,
                          minReps = 7,
                          sortColumns = F,
                          column_sort_factors = NA,
                          customAnno = NULL,
                          annoJoinCol = NA,
                          baseMeanCutoff = 0,
                          random_seed = 12345,
                          custom_order = NULL,
                          tidy_go = FALSE,
                          return_fold_enrichment = FALSE,
                          ...)
{
  library(pheatmap)
  library(RColorBrewer)
  library(topGO)
  library(GO.db)
  library(DESeq2)
  set.seed(random_seed)
  
  counts_mat = construct_goi_matrix(dge = dge,
                                    qval_cutoff = qval_cutoff,
                                    genes_of_interest = genes_of_interest,
                                    design = design,
                                    cores = cores,
                                    minReps = minReps)
  
  set.seed(random_seed)
  kmeans_results = as.data.frame(kmeans(x = counts_mat,
                          centers = k, 
                          iter.max = 1000, 
                          nstart = 25)$cluster)
  kmeans_results = rownames_to_column(kmeans_results,
                                      var = "gene")
  colnames(kmeans_results)[2] = "cluster"
  #if user specifies a custom order of clusters (top to bottom)
  if(!is.null(custom_order))
  {
    cluster_conv = data.frame(new_cluster = c(1:max(kmeans_results$cluster)),
                              cluster = custom_order)
    kmeans_results = kmeans_results %>%
      dplyr::left_join(., cluster_conv) %>% 
      dplyr::arrange(new_cluster)
    colnames(kmeans_results) = c("gene", "old_cluster", "cluster") #fit into following code
  }
  
  #order genes by cluster assignment
  kmeans_results = kmeans_results[order(kmeans_results$cluster), ]
  counts_mat = counts_mat[kmeans_results$gene, ]
  
  #if requested, sort columns by any number of factors
  md = as.data.frame(colData(dge))
  if(sortColumns)
  {
    columns_sorted = md %>% 
      dplyr::arrange_(.dots = column_sort_factors)
    counts_mat = counts_mat[, rownames(columns_sorted)]
  }
  
  #generate gaps for each cluster
  gaps = c()
  for(i in 2:nrow(kmeans_results))
  {
    prev = kmeans_results$cluster[i - 1]
    cur = kmeans_results$cluster[i]
    if(cur != prev)
    {
      gaps = append(gaps, (i - 1))
    }
  }
  
  # extract metadata for legend (if necessary)
  if(!is.null(legend_factors))
  {
    md = md[, legend_factors, drop = F]
  }
  
  #add GO terms, as necessary
  if(display_go_terms || return_go_terms)
  {
    library(go_annotations, character.only = T) #load GO package using variable
    #define universe as all detected genes in dataset
    all_counts = counts(dge, normalized = T)
    universe = rownames(all_counts[rowSums(all_counts) > 0, ])
    fisherTest = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
    
    #determine gene ID type automatically
    #assume some transgenes will not fit the regular expression
    if(sum(grepl("^ENS|WBcel", rownames(dge))) > 100)
    {
      message("Using ID type: ensembl_gene_id")
      id_type = "ensembl"
      #entrezgenes are all numeric
    } else if(sum(grepl("^\\d+$", rownames(dge))) > 100)
    {
      message("Using ID type: entrezgene_id")
      id_type = "entrez"
    } else
    {
      message("Using ID type: external_gene_name")
      id_type = "genename"
    }
    
    clusters = unique(as.character(kmeans_results$cluster))
    cluster_GO = mclapply(clusters, function(x){
      cluster_genes = kmeans_results[kmeans_results$cluster == x, "gene"]
      selection = as.numeric(universe %in% cluster_genes)
      names(selection) = universe
      go_data = new("topGOdata", 
                    ontology = go_ontology, 
                    allGenes = selection,
                    geneSel = function(x){
                      return(x == 1)},
                    annot = annFUN.org, 
                    mapping = go_annotations, 
                    ID = id_type)
      
      #run Fisher test
      test_results = topGO::getSigGroups(go_data, fisherTest)
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
        score = rownames_to_column(score, var = "go_id")
      
        #adjust p-values and take significant
        score$padj = p.adjust(score$pval, method = "fdr")
        score = subset(score, padj < 0.05)
      }
      #in case of no significant go terms, return NULL
      if(nrow(score) == 0)
      {
        return(NULL)
      } else
      {
        score$description = NA
        for(i in 1:nrow(score))
        {
          score$description[i] = GOTERM[[score$go_id[i]]]@Term
        }
        
        #add descriptions
        score$full_go = paste(score$go_id, score$description)
        return(score)
       }},
      mc.preschedule = T, 
      mc.cores = getOption("mc.cores", cores))
    
    if(display_go_terms)
    {
      #annotate clusters using empty gene name slots
      cluster_annos = rep("", nrow(counts_mat))
      for(i in 1:length(cluster_GO))
      {
        #skip empties
        if(is.null(cluster_GO[[i]]))
        {
             next
        }
        if(i == 1)
        {
          start = 1
        } else
        {
          start = gaps[i - 1] + 1
        }
        
        #sort by p-value
        cluster_GO[[i]] = cluster_GO[[i]][order(cluster_GO[[i]]$padj), ]
        
        #add spacing to make it readable
        if(nrow(cluster_GO[[i]]) >= max_go_terms)
        {
          n_go = max_go_terms
        } else
        {
          n_go = nrow(cluster_GO[[i]])
        }
        
        #space based on size of cluster (fill just top half)
        skip = max((nrow(cluster_GO[[i]]) / max_go_terms / 2), 75) #minimum of 75 for this font
        for(j in 1:n_go)
        {
          cur_index = start + skip * (j - 1)
          #skip term if we've filled the cluster space already
          cluster_end = ifelse(i == length(cluster_GO),
                               yes = nrow(counts_mat),
                               no = gaps[i] - 1) #last cluster doesn't have a following start
          if(cur_index <= cluster_end)
          {
            cluster_annos[cur_index] = cluster_GO[[i]]$full_go[j]
          }
        }
      }
    }
  }
  if(!is.null(customAnno) && display_go_terms == FALSE)
  {
    #pare down to genes in the matrix
    customAnno = tibble::column_to_rownames(customAnno, var = annoJoinCol)
    customAnno = customAnno[rownames(counts_mat), ]
    cluster_annos = customAnno
  } else if(is.null(customAnno) && display_go_terms == FALSE)
  {
    cluster_annos = NULL
  }
  plot = pheatmap::pheatmap(counts_mat, 
                  cluster_rows = F,
                  cluster_cols = cluster_columns,
                  clustering_method = "ward.D2",
                  gaps_row = gaps,
                  show_colnames = F,
                  annotation_col = md,
                  labels_row = cluster_annos,
                  fontsize_row = label_fontsize,
                  annotation_names_col = F,
                  color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                          name = "RdBu")))(100),
                  ...)
  
  output = list("plot" = plot, "genes" = NULL, "GO" = NULL)
  if(return_genes)
  {
    library(biomaRt)
    mart = useMart("ensembl", ensembl_db)
    conv = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                 mart = mart)
    kmeans_results = dplyr::left_join(kmeans_results,
                               conv,
                           by = c("gene" = "ensembl_gene_id")) %>% 
      dplyr::arrange(cluster)
    output$genes = kmeans_results
  }
  if(return_go_terms)
  {
    if(tidy_go) #join up into data frame
    {
      for(i in 1:length(cluster_GO))
      {
        if(!is.null(cluster_GO[[i]]))
        {
          cluster_GO[[i]]$cluster = i
        }
      }
      #catch error when there are no terms at all
      non_null = sum(vapply(cluster_GO, function(x){ !is.null(x) }, FUN.VALUE = 1))
      if(non_null > 0)
      {
        cluster_GO = bind_rows(cluster_GO) %>% 
          arrange(cluster, padj)
      } else
      {
        message("No significant go terms to return")
      }
    }
    non_null = sum(vapply(cluster_GO, function(x){ !is.null(x) }, FUN.VALUE = 1))
    if(non_null == 0)
    {
      message("No significant go terms to return")
    }
    output$GO = cluster_GO
  }
  
  return(output)
}
  
  