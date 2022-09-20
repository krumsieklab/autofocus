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
#'
#' @return Returns the densities and number of significant descendants of each node
#'
#' @noRd
get_sig_child_density <- function(R, phenotype, confounders){
  sig_kids<-lapply(1:dendextend::nleaves(R$HCL), function(i) {
    data <- data.frame(pheno=R$samples[[phenotype]], molecule = R$data[,i], R$samples[,confounders])
    reg <- lm(pheno ~ ., data = data)
    summary(reg)$coefficients["molecule", "Pr(>|t|)"]

  }) %>% unlist %>%  p.adjust(method="bonferroni")

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
