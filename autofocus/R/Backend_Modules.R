#' Get dendogram
#'
#' Hierarchically clusters data in R.
#'
#' @param R R struct.
#' @param method Hierarchical clustering method. If NA (default), get optimal choice of all clustering methods.
#' @param cores Number of cores to be used in this step.
#'
#' @return hclust object of dendrogram
#'
#' @noRd
get_dendro <- function(R,
                       method = NA,
                       cores = cores
){

  # Calculate distance metric
  dist_mat <- abs_cor_dist(R)

  # Hclust using defined method if not NA
  if (!(is.na(method))){
    return(hclust(dist_mat, method = method))
  }

  # Hclust using search over all clustering methods
  else {
    return(get_cluster_method(dist_mat, cores = cores))
  }
}

#' Distant metric
#'
#' Calculates distance matrix of data based on correlation value
#'
#' @param R R struct.
#'
#' @return Distance matrix made up of 1 - the absolution value of the correlation matrix
#'
#' @noRd
abs_cor_dist <- function(R){
  as.dist((1-abs(R$C))) # abs correlation
}

#' Get indices in dendogram
#'
#' Recursively finds the indices of each node in the dendrograms "merge" object.
#'
#' @param R R struct.
#' @param number Number of the split in the dendrogram.
#' @param indices List to build the indices.
#'
#' @return the indices of nodes in R's hclust object
#'
#' @noRd
get_dend_indices <- function(R, number, indices){
  me <- R$HCL$merge

  #Leaf case
  if(number <0){
    indices <- c(indices, abs(number) + dim(me)[1])
  }

  else{
    left_indices <- get_dend_indices(R, me[number,][1], indices)
    right_indices <-get_dend_indices(R, me[number,][2], indices)
    indices <- c(indices, number, left_indices, right_indices)
  }
  indices
}

#' Get Cluster members
#'
#' Recursively collects the indices of leaves that are descendents of a module.
#'
#' @param R R struct.
#' @param i Index of module we are finding the members of.
#' @param cl Vector to collect node indices.
#'
#' @return vector of indices of leaves in a module
#'
#' @noRd
get_members <- function(R,
                        i,
                        cl = c()
){
  if (i<0){
    return(i)
  }
  x <- R$HCL$merge[i,]
  if(sum(x<0)==2) return(c(cl,-x))
  if(x[1]<0 & x[2]>0) return(get_members(R, x[2],c(cl,-x[1])))
  if(x[2]<0 & x[1]>0) return(get_members(R,x[1],c(cl,-x[2])))
  return(c(cl,get_members(R,x[1]),get_members(R,x[2])))
}

#' Get parent of a cluster
#'
#' Recursively collects the indices of
#' leaves that are descendents of a module
#'
#' @param R R struct
#' @param i index of module we are finding the members of
#' @param cl vector to collect node indices
#'
#' @return vector of indices of leaves in a module
#'
#' @noRd
get_parent <- function(
    R,
    i
){
  if (i > dim(R$HCL$merge)[1]){
    ind <- dim(R$HCL$merge)[1] - i
    return(max(which(R$HCL$merge[,2]  == ind), which(R$HCL$merge[,1]  == ind)))
  }
  max(which(R$HCL$merge[,2]  == i), which(R$HCL$merge[,1]  == i))
}

#' Get proportion/density of significant children
#'
#' Finds the proportion of sigificant leaves at each
#' parent node of the dendrogram
#'
#' @param R R struct
#' @param phenotype disease pheotype
#' @param confounders vectors of confounder names
#' @param method Method to use for p-value adjustment. Default: "fdr".
#'
#' @return Returns the densities and number of significant descendants of each node
#'
#' @noRd
get_sig_child_density <- function(R, phenotype, confounders, method="fdr"){
  sig_kids<-lapply(1:dendextend::nleaves(R$HCL), function(i) {
    data <- data.frame(pheno=R$samples[[phenotype]], molecule = R$data[,i], R$samples[,confounders])
    reg <- lm(pheno ~ ., data = data)
    summary(reg)$coefficients["molecule", "Pr(>|t|)"]

  }) %>% unlist %>%  p.adjust(method=method)

  sig_kids<-1*(sig_kids<0.05)

  densities <- lapply(1:dim(R$HCL$merge)[1], function(i) {
    members <-  R$clusts[i][[1]]
    sum(sig_kids[members])/length(members)
  }) %>% unlist()
  c(densities, sig_kids)
}

#' Get direct edges between disease phenotype and module members
#'
#' Uses mixed graphical models to create a graph
#' of direct edges between module members, confounders,
#' and disease phenotype
#'
#' @param i index of parent node of cluster
#' @param R R struct
#' @param phenotype disease pheotype
#' @param confounders vectors of confounder names
#'
#' @return Returns a graph object
#'
#' @noRd
get_edges_linear <- function(i, R, phenotype, confounders){

  if(length(phenotype)==1){
    if (!(R$clust_info$mgm_capable[i])){
      return (NA)
    }

    else{
      members <-  R$clusts[i][[1]]

      # Get variable type of phenotype and confounders
      # Either categorical ("c") or gaussian ("g")
      # If categorical, get number of categories

      #phenotype_type <- get_variable_type(R$samples[[phenotype]])
      #phenotype_cat <- get_num_categories(R$samples[[phenotype]])

      confounder_data <- R$samples[,confounders]


      #Combine all data with confounders
      full_data<-data.frame(phenotype = R$samples[[phenotype]],
                            R$data[,members],
                            confounder_data)

      full_data <-as.matrix(as.data.frame(lapply(full_data, as.numeric)))
      # Remove samples with NAs
      no_nas = apply(is.na(full_data),1,sum) == 0
      data_no_na <- full_data[no_nas,]

      types <- apply(data_no_na, 2, get_variable_type) %>% as.vector()
      levels <- apply(data_no_na, 2, get_num_categories) %>% as.vector()

      if('c' %in% types) {
        ind_cat <- which(types == 'c')
        to_remove<-c()
        for(i in ind_cat) {
          l_frqu <- table(data_no_na[,i])
          cats_to_remove <- names(l_frqu)[which(l_frqu<=1)]
          if(length(cats_to_remove)>0){
            levels[i]<-levels[i] - length(cats_to_remove)
            to_remove<-c(to_remove, which(data_no_na[,i]%in% cats_to_remove))
          }
        } # this does not catch the case where one category is not present at all; but this is catched by comparing specified levels and real levels
      }
      if(length(to_remove)>0){
        data_no_na <- data_no_na[-to_remove,]
      }
      # Find direct connections with phenotype in mgm
      fit_mgm <- mgm::mgm(data_no_na, type=types, level=levels, lambdaSel="EBIC")
      adj_mat <- 1* (fit_mgm$pairwise$wadj!=0)

      # Create graph
      colnames(adj_mat)<- c(phenotype, colnames(R$data)[members], confounders)
      G<-igraph::graph_from_adjacency_matrix(adj_mat, add.rownames = T)
      igraph::V(G)$NodeType <- c("phenotype", rep("analyte",length(members)), rep("confounder", length(confounders)))
      igraph::V(G)$Driver <- igraph::V(G)$name %in% igraph::neighbors(G,phenotype)$name[igraph::neighbors(G,phenotype)$NodeType =="analyte"]
      igraph::V(G)$color <- c("blue", rep("red",length(members)), rep("grey", length(confounders)))
      G
    }
  }else if(length(phenotype)==2){
    if (!(R$clust_info$pheno1_mgm_capable[i]) & !(R$clust_info$pheno2_mgm_capable[i])){
      return (NA)
    }

    else{
      members <-  R$clusts[i][[1]]

      # Get variable type of phenotype and confounders
      # Either categorical ("c") or gaussian ("g")
      # If categorical, get number of categories

      #phenotype_type <- get_variable_type(R$samples[[phenotype]])
      #phenotype_cat <- get_num_categories(R$samples[[phenotype]])

      confounder_data <- R$samples[,confounders]


      #Combine all data with confounders
      full_data<-data.frame(R$data[,members],
                            confounder_data)
      full_data$phenotype1 = R$samples[[phenotype[1]]]
      full_data$phenotype2 = R$samples[[phenotype[2]]]

      pheno_names <- phenotype

      full_data <-as.matrix(as.data.frame(lapply(full_data, as.numeric)))
      # Remove samples with NAs
      no_nas = apply(is.na(full_data),1,sum) == 0
      data_no_na <- full_data[no_nas,]

      types <- apply(data_no_na, 2, get_variable_type) %>% as.vector()
      levels <- apply(data_no_na, 2, get_num_categories) %>% as.vector()

      if('c' %in% types) {
        ind_cat <- which(types == 'c')
        to_remove<-c()
        for(i in ind_cat) {
          l_frqu <- table(data_no_na[,i])
          cats_to_remove <- names(l_frqu)[which(l_frqu<=1)]
          if(length(cats_to_remove)>0){
            levels[i]<-levels[i] - length(cats_to_remove)
            to_remove<-c(to_remove, which(data_no_na[,i]%in% cats_to_remove))
          }
        } # this does not catch the case where one category is not present at all; but this is catched by comparing specified levels and real levels
      }
      if(length(to_remove)>0){
        data_no_na <- data_no_na[-to_remove,]
      }
      # Find direct connections with phenotype in mgm
      fit_mgm <- mgm::mgm(data_no_na, type=types, level=levels, lambdaSel="EBIC")
      adj_mat <- 1* (fit_mgm$pairwise$wadj!=0)

      # Create graph
      colnames(adj_mat)<- c(colnames(R$data)[members], confounders, pheno_names)
      G<-igraph::graph_from_adjacency_matrix(adj_mat, add.rownames = T)
      igraph::V(G)$NodeType <- c(rep("analyte",length(members)), rep("confounder", length(confounders)), rep("phenotype", length(pheno_names)))
      igraph::V(G)$Driver <- igraph::V(G)$name %in% igraph::neighbors(G,pheno_names)$name[igraph::neighbors(G,pheno_names)$NodeType =="analyte"]
      igraph::V(G)$color <- c(rep("red",length(members)), rep("grey", length(confounders)), rep("blue",length(pheno_names)))
      G
    }
  }else{
    stop("autofocus does not support more than two phenotypes.")
  }
}

#' Determine whether a variable is discrete or continuous
#'
#' Determines whether a variable is discrete or continuous based
#' on the ratio of the number of unique entries to number of entries
#'
#' @param variable_vec The vector to determine type
#' @param cutoff The cutoff of the unique entries to length ratio
#'
#' @return Returns "c" if discrete categorial or "g" if continuous or gaussian
#'
#' @noRd
get_variable_type <- function(variable_vec, cutoff= 0.05){
  ifelse((length(unique(variable_vec))/length(variable_vec)) <= cutoff, "c", "g")
}

#' Determine number of categories in a variable
#'
#' If a variable is determinned to be categorical, the number of categories
#' Else returns 1
#'
#' @param variable_vec The vector to determine category count
#' @param cutoff The cutoff of the unique entries to length ratio
#'
#' @return Returns number of categories or 1 if continuous or Gaussian
#'
#' @noRd
get_num_categories <- function(variable_vec, cutoff=0.05){
  ifelse((length(unique(variable_vec))/length(variable_vec)) <= 0.05, length(unique(variable_vec)), 1)

}

#' Use WGCNA to get adjacency matrix
#'
#' Helper function called by get_denro that creates an adjacency matrix. The matrix is created
#' using the saved Topological Overlap Matrices (TOMs) saved by the WGCNA function
#' blockwiseModules.
#'
#' @param input_mat R Data matrix with samples in rows and features in columns.
#' @param cores Number of cores to be used in this step.
#'
#' @return NULL (save TOM file to working directory)
#'
#' @import WGCNA
#'
#' @noRd
run_wgcna <- function(input_mat,
                      cores
){
  allowWGCNAThreads(nThreads = cores)

  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

  # Call the network topology analysis function
  sft = WGCNA::pickSoftThreshold(
    input_mat,
    powerVector = powers,
    verbose = 5
  )

  # Select the power with the best fit
  signed_r2 = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
  select_power = which(signed_r2 == max(signed_r2))

  # module detection (use WGCNA cor implementation for this step)
  cor <- WGCNA::cor
  netwk <- WGCNA::blockwiseModules(input_mat,                # <= input here
                            # == Adjacency Function ==
                            power = select_power, # <= power here
                            blocks = rep(1, ncol(input_mat)),

                            # == Tree and Block Options ==
                            minModuleSize = 2,
                            saveTOMs = T,

                            # == Output Options
                            numericLabels = T,
                            verbose = 3)
  cor <- stats::cor

  # the above WGCNA function saves the TOM file to the local directory
  # so nothing is returned from this function
  return(NULL)
}

#' Get WGCNA TOM-Based Distance Metric
#'
#' Calculate the WGCNA Distance Metric equal to one minus the topological overlap matrix. Calls
#' helper function run_wgcna. See \code{\link[WGCNA]{blockwiseModules}} for details.
#'
#' @param input_matrix R Data matrix with samples in rows and features in columns.
#' @param cores Number of cores to be used in this step.
#'
#' @return TOM-Based Distance Metric
#'
#' @noRd
get_wgcna_dist_metric <- function(input_matrix, cores){

  # calculate TOM
  run_wgcna(input_matrix, cores)

  # load TOM file and calculate distance metric
  tom_file <- "blockwiseTOM-block.1.RData"
  load(tom_file)
  tom = as.matrix(TOM)
  rm(TOM)
  WGCNA::collectGarbage()
  dissTom = 1 - tom

  # clean up and remove temporary TOM file
  rm(tom)
  WGCNA::collectGarbage()
  file.remove(tom_file)

  return(dissTom)

}

#' Get summary cluster information at an input threshold
#'
#' This is a helper function for the threshold_analysis function. It calculates
#' the metrics for threshold analysis for an individual threshold
#'
#' @param R R_struct output from initialize_R
#' @param threshold Which threshold to calculate metric for
#' @param density_col Column name in R struct that stores densities
#' @param metric Which metric to calculate
#'
#' @return The input metric for the clusters returned by R at input threshold
#'
#' @noRd
get_cluster_info <- function(
    R,
    threshold,
    metric = "prop_non_sig_in_clusts"){
  # Get tree
  hc <- R$HCL

  # Single Phenotype
  # Get all peaks meeting threshold
  potential_peaks <- identify_peaks(R, dim(hc$merge)[1], threshold, c(), "densities")
  # Filter out piggy backers, replaced with their children
  filtered_peaks <- filter_peaks(R, potential_peaks, threshold, "densities")

  # While there are still new children to check
  while(!(all(potential_peaks==filtered_peaks))){

    # Get the unchecked children
    new_potential_peaks <- unlist(lapply(filtered_peaks,function(i)
      identify_peaks(R, i, threshold, c(), "densities")))
    filtered_peaks <- filter_peaks(R, new_potential_peaks, threshold, "densities")
    potential_peaks <- new_potential_peaks
  }
  nodes_in_clusters <- R$clust_info %>%
    mutate(num_sig_in_clust = Size*densities) %>%
    mutate(num_non_sig_in_clust = Size - (Size*densities)) %>%
    filter(ClusterID %in% filtered_peaks)

  # Proportion of nodes in returned clusters that aren't significant
  if(metric == "prop_non_sig_in_clusts"){
    return(sum(nodes_in_clusters$num_non_sig_in_clust)/sum(nodes_in_clusters$Size))
  }
  # Number of singletons
  else if(metric =="singletons"){
    return(sum(R$clust_info$Size==1 & R$clust_info$densities==1) - sum(nodes_in_clusters$num_sig_in_clust))
  }
  # Number of clusters
  else if(metric == "num_clusts"){
    return(length(filtered_peaks))
  }
  #  Cluster size range
  else if(metric == "clust_size_range"){
    return(nodes_in_clusters %>% .$Size %>% range %>% diff)
  }
  # Average Enrichment of clusters
  else if(metric == "avg_enrichment"){
    return(mean(nodes_in_clusters$densities))
  }
  else{ # Height range
    nodes_in_clusters %>% .$Coord_Y %>% range %>% diff
  }
}
