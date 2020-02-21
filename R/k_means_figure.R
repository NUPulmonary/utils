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

#for extracting counts for genes of interest from a DESeqDataSet
construct_goi_matrix = function(dge,
                                qval_cutoff = 0.05,
                                genes_of_interest = NULL,
                                design = NULL,
                                cores = 1,
                                fitType = "local")
{
  register(MulticoreParam(cores))
  
  counts_mat = counts(dge, normalized = T)
  if(is.na(qval_cutoff) || (is.null(design) && is.null(genes_of_interest))) #all cases to include all genes
  {
    qval_cutoff = 1
  }
  
  #determine subset of genes to use
  if(is.null(genes_of_interest) && qval_cutoff < 1)
  {
    #perform anova-like test to identify genes which vary significantly across factor(s) of interest
    design = as.formula(design)
    design(dge) = design
    deseq_results = DESeq(dge, 
                          test = "LRT", 
                          reduced = ~ 1,
                          parallel = T,
                          fitType = fitType)
    deseq_results = as.data.frame(results(deseq_results))
    genes_of_interest = rownames(subset(deseq_results, padj < qval_cutoff))
  } else
  {
    genes_of_interest = rownames(counts_mat)
  }
  
  #subset counts down to chosen genes
  counts_mat = counts_mat[genes_of_interest, ]
  
  #z-score matrix (necessary for k means)
  counts_mat = t(scale(t(counts_mat))) #transpose issue is annoying
  
  return(counts_mat)
}
  

#for heuristically determining an ideal k (run this first)
k_elbow = function(dge,
                   qval_cutoff = 0.05,
                   genes_of_interest = NULL,
                   design = NA,
                   cores = 1,
                   max_k = 50)
{
  require(ggplot2)
  require(tibble)
  
  counts_mat = construct_goi_matrix(dge = dge,
                                    qval_cutoff = qval_cutoff,
                                    genes_of_interest = genes_of_interest,
                                    design = design,
                                    cores = cores)
  
  #now run kmeans for all values of k, and find sums of squared differences within each cluster for each k
  sums_of_squares = mclapply(1:max_k, function(k){
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
                          max_go_terms = 5,
                          colnames = F,
                          legend_factors = NULL,
                          go_annotations = "org.Mm.eg.db",
                          mart_name = "mmusculus_gene_ensembl")
{
  require(pheatmap)
  require(RColorBrewer)
  require(topGO)
  require(GO.db)
  
  counts_mat = construct_goi_matrix(dge = dge,
                                    qval_cutoff = qval_cutoff,
                                    genes_of_interest = genes_of_interest,
                                    design = design,
                                    cores = cores)
  
  kmeans_results = as.data.frame(kmeans(x = counts_mat,
                          centers = k, 
                          iter.max = 1000)$cluster,
                          nstart = 25)
  kmeans_results = rownames_to_column(kmeans_results,
                                      var = "gene")
  colnames(kmeans_results)[2] = "cluster"
  
  #order genes by cluster assignment
  kmeans_results = kmeans_results[order(kmeans_results$cluster), ]
  counts_mat = counts_mat[kmeans_results$gene, ]
  
  #generate breaks for each cluster
  breaks = c()
  for(i in 2:nrow(kmeans_results))
  {
    prev = kmeans_results$cluster[i - 1]
    cur = kmeans_results$cluster[i]
    if(cur != prev)
    {
      breaks = append(breaks, i)
    }
  }
  
  # extract metadata for legend (if necessary)
  if(!is.null(legend_factors))
  {
    md = as.data.frame(colData(dge))
    md = md[, legend_factors]
  }
  
  #add GO terms, as necessary
  if(display_go_terms)
  {
    #define universe as all detected genes in dataset
    all_counts = counts(dge, normalized = T)
    universe = rownames(all_counts[rowSums(all_counts) > 0, ])
    fisherTest = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
    
    clusters = unique(as.character(kmeans_results$cluster))
    cluster_GO = lapply(clusters, function(x){
      cluster_genes = kmeans_results[kmeans_results$cluster == x, "gene"]
      selection = as.numeric(universe %in% cluster_genes)
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
      score = rownames_to_column(score, var = "go_id")
      
      #adjust p-values and take significant
      score$padj = p.adjust(score$pval, method = "fdr")
      score = subset(score, padj < 0.05)
      score$description = NA
      for(i in 1:nrow(score))
      {
        score$description[i] = GOTERM[[score$go_id[i]]]@Term
      }
      
      #add descriptions
      score$full_go = paste(score$go_id, score$description)
      return(score)})
      
    #annotate clusters using empty gene name slots
    cluster_annos = rep("", nrow(counts_mat))
    for(i in 1:length(cluster_GO))
    {
      if(i == 1)
      {
        start = 1
      } else
      {
        start = breaks[i - 1] + 1
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
      for(j in 1:n_go)
      {
        cur_index = start + 100 * (j - 1)
        cluster_annos[cur_index] = cluster_GO[[i]]$full_go[j]
      }
    }
  }
  plot = pheatmap(counts_mat, 
                  cluster_rows = F,
                  cluster_cols = T,
                  clustering_method = "ward.D2",
                  gaps_row = breaks,
                  show_colnames = F,
                  annotation_col = md,
                  labels_row = cluster_annos,
                  fontsize_row = 6,
                  annotation_names_col = F,
                  color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                          name = "RdBu")))(100))
  return(plot)
}
  
  