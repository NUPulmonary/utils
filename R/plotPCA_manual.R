# copy of DESeq2 plotPCA, but returns the actual PCA object
# add-ins
# now returns percent variance for each PC requested
# merge_metadata: joins output dataframe with colData. Defaults to FALSE for compatibility reasons.
# return_loadings: return gene loadings data frame. Defaults to FALSE for compatibility reasons.
# mart_name: prefix for biomart dataset name. Defaults to "mmusculus". Ignored if custom annotation is supplied or return_loadings is FALSE.
# custom_annotation: a dataframe with "ensembl_gene_id" and "external_gene_name" as columns. Take precedence over mart name if supplied.
plotPCA_manual = function(object, 
                          intgroup="condition", #essentially ignored if merge_metadata = TRUE
                          ntop=500, 
                          pcs = 2,
                          merge_metadata = FALSE,
                          return_loadings = FALSE,
                          mart_name = "mmusculus",
                          custom_annotation = NULL)
{
  require(tidyverse)
  
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assemble the data for the plot
  d <- data.frame(group=group, intgroup.df, name=colnames(object))
  tmp = pca$x[, 1:pcs]
  colnames(tmp) = paste0("PC", 1:pcs)
  d = cbind(tmp, d)
  attr(d, "percentVar") <- percentVar[1:pcs]
  
  if(merge_metadata == TRUE)
  {
    md = colData(object) %>% 
      as.data.frame() %>% 
      rownames_to_column("sample")
    d = d %>% 
      dplyr::select(-group) %>% 
      dplyr::rename(sample = name) %>% 
      left_join(., md) %>% 
      column_to_rownames("sample")
  }
      
  
  plt = ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed()
  
  #if requested, return loadings in human-readable format
  if(return_loadings == TRUE)
  {
    if(is.null(custom_annotation))
    {
      require(biomaRt)
      mart = useMart("ensembl", paste0(mart_name, "_gene_ensembl"))
      conv = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   mart = mart)
    } else
    {
      conv = custom_annotation
    }
    
    loadings = pca$rotation %>% 
      as.data.frame() %>% 
      rownames_to_column("ensembl_gene_id") %>% 
      left_join(., conv) %>% 
      dplyr::relocate(ensembl_gene_id, external_gene_name) %>% 
      arrange(desc(PC1))
  } else
  {
    loadings = NULL
  }
  
  #output percent variance for each PC as a tidy data frame
  percent_var = data.frame(component = 1:pcs, percent_var = round(percentVar[1:pcs] * 100))
  
  return(list(pca = pca, data = d, plot = plt, loadings = loadings, percent_var = percent_var))
}
