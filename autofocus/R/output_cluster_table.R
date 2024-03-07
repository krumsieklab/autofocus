#' Output AutoFocus clusters at a specified threshold in table format
#'
#' Outputs a dataframe containing the information on nodes found in AutoFocus clusters of size greater than 1
#'
#' @param R R_struct object output by initialize_R
#' @param threshold Enrichment threshold at which to calculate cluster enrichment
#' @param phenotype Which phenotype associated clusters to look at. Default "both" to see all phenotypes' clusters
#' @param name Column name in R$annos specifying the node name
#' @param annotations_to_show Vector of colnames of R$annos to return with nodes
#'
#' @return Dataframe with the nodes in clusters, their p-value association with phenotypes of interest, cluster, and annotations_to_show
#'
#'
#' @author AS
#'
#'
#' @export
output_cluster_table <- function(R, threshold, phenotype = "both", name, annotations_to_show){
  peaks =  peak_finder_wrapper(R, threshold = threshold)
  clusts_to_show <- R$clust_info[(peaks!="#666666"&R$clust_info$Size>1 ),] # Get the clusters that are peaks and not leaves
  nodes_to_show <- R$clust_info

  # this subselects the annos for the clusters of interest (from R$annos)
  # and adds a column that just assigns the cluster number to each molecule involved
  lapply(1:nrow(clusts_to_show), function(i){
    annos <- R$annos[R$clusts[[clusts_to_show$ClusterID[i]]],] # Get annotations for molecules in the peak clusters
    annos$cluster =rep(clusts_to_show$ClusterID[i], nrow(annos)) # Make a new column that just reps the cluster id for each molecule
    data.frame(annos[,c(name, annotations_to_show, "p_value","cluster")])

  }) %>% bind_rows()

}
