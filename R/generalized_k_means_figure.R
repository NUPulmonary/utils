#' Generalized function for generating a clustermap
#' with rows split into clusters using k-means
#' 
#' @param mat matrix of values where rows are features and columns are samples
#' @param md dataframe where row names are samples and columns represent associated metadata
#' @param k chosen value of k for k-means
#' @return a list object containing 'hm': the heatmap and 'assignments': cluster assignments for each feature
#' @export

#' Function to estimate ideal k using the within-cluster sum-of-squares method
#' 
#' @param mat matrix of values where rows are features and columns are samples
#' @param max_k maximum value of k to use for plotting; defaults to 50
#' @param cores number of cores to use for simultaneous estimation of values of k; defaults to 1
#' @param random_seed random seed for kmeans function
k_elbow = function(mat,
                   random_seed = 12345,
                   cores = 8,
                   max_k = 50)
{
  require(BiocParallel)
  register(MulticoreParam(cores))
  
  #now run kmeans for all values of k, and find sums of squared differences within each cluster for each k
  sums_of_squares = mclapply(1:max_k, function(k){
    set.seed(random_seed)
    kmeans_result = kmeans(x = na.omit(mat), 
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