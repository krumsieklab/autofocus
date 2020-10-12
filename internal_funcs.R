### internal functions

#source("~/Documents/JanLab/git/autofocus/Backend_Modules.R")

#### Initialize ####
#' initialize
#' 
#' Takes as input a data matrix, annotations on the samples and annotations on the molecules
#' and creates a list structure, R with the hierarchical clustering structure and cluster informations
#'

#' @param data.matrix Matrix with samples as rows and biomolecules as columns
#' @param sample.data Matrix of annotation data of samples
#' @param mol.data Matrix of annotation data of biomolecules
#' 
#' @return R struct with the data matrix, correlation matrix, 
#' sample and molecular annotations, hierarchical structure, cluster membership
#' node order within the hierarchical structure and a color vector for the platforms

initialize <- function(
  data.matrix, # Matrix of raw data (samples * nodes)
  sample.data, # Matrix of sample annotations
  mol.data, # Matrix of node annotations
  cores
){
  
  R <- list(data = data.matrix, samples = sample.data, annos = mol.data)
  R$C <- cor(data.matrix)
  
  R$HCL <- R %>% get_dendro(method = "average", cores)
  R$order <- get_dend_indices(R, dim(R$HCL$merge)[1], c())
  
  R$platforms = unique(R$annos$platform)
  R$clusts <- lapply(1:dim(R$HCL$merge)[1], function(x) get_members(R, x))
  
  R
}

#### Color Nodes ####
#' get_node_color
#' 
#' Gets the color of the nodes within the hierarchical tree based on 
#' BIC and p-value
#'
#' @param R R struct
#' @param i index of module we are coloring
#' @param signif Is the module i significant?
#' 
#' @return color of node i 
#
get_node_color <- function(
  R,
  i,
  signif,
  node_color,
  node_color_light
){
  hc <- R$HCL
  internal_nodes <- dim(hc$merge)[1]
  if (i %in% signif[!is.na(signif)]){
    if (i > (internal_nodes)) return(node_color)
    else {
      
      children <- hc$merge[i,]
      child1 <- if (children[1]<0) (internal_nodes + abs(children[1])) else children[1]
      child2 <- if (children[2]<0) (internal_nodes + abs(children[2])) else children[2]
      
      ## Children aren't significant, combination is
      if (R$pvals[child1] > 0.05 & R$pvals[child2] > 0.05) return(node_color)
      
      ## Children are significant
      if(R$pvals[child1] < 0.05 & R$pvals[child2] < 0.05){
        
        if (R$BIC[[i]] < R$BIC[[child1]] & R$BIC[[i]]<R$BIC[[child2]]) return(node_color) else return(node_color_light)
      }
      
      ## One child is significant
      else {
        sig_child <- if (R$pvals[child1] < R$pvals[child2]) child1 else child2
        
        if(R$BIC[[i]] < R$BIC[[sig_child]]) return(node_color) else return(node_color_light)
      }
    }
  }
  else{
    return("black")
  }
}


#### Get Cluster members ####
#' get_members
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
get_members <- function(
  R,
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





#### Convert dendrogram to adjacency matrix ####
#' dend_to_adj_mat
#' 
#' Convert dendrogram to adjacency matrix for visualization
#'
#' @param hclust_obj the hclust object to convert
#' 
#' @return an adjacency matrix of the dendrogram where each parent has an edge with its children
#' 

dend_to_adj_mat <- function(
  hclust_obj
){
  
  adj_mat <- matrix(0L, nrow = (2*dim(hclust_obj$merge)[1])+1, ncol = (2*dim(hclust_obj$merge)[1])+1)
  for(i in rev(1:dim(hclust_obj$merge)[1])){
    parent <- i
    right_index <- hclust_obj$merge[i,1]
    if (right_index < 0){
      right_index <- abs(right_index) + dim(hclust_obj$merge)[1]
    }
    left_index <- hclust_obj$merge[i,2]
    if (left_index < 0){
      left_index <- abs(left_index) + dim(hclust_obj$merge)[1]
    }
    adj_mat[i, right_index] <- 1
    adj_mat[i, left_index] <- 1
  }
  adj_mat
}

#### Get indices in dendrogram ####
#' get_dend_indices
#' 
#' Recursively finds the indices of each node in the dendrograms "merge" object 
#'
#' @param R R struct
#' @param number number of the split in the dendrogram
#' @param indices list to build the indices
#' 
#' @return the indices of nodes in R's hclust object
#' 
get_dend_indices <- function(
  R,
  number, 
  indices
){
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

#### Output Structure ####
#' make_R
#' 
#' Wrapper function taking in preprocessed data, outputting significantly labelled dendrogram
#'
#' @param SE SummarizedExperiment object with preprocessed data
#' @param phenotype name of the phenotype we are testing for association
#' @param confounders vector of confounder names
#' @param save_file do you want to save the R struct?
#' @param filename name of file if save_file is true
#' @param cores number of cores to use to run algorithm
#' 
#' @return R struct with significance data to be passed into shiny app
#' 

make_R <- function(
  SE,
  phenotype,
  confounders,
  node_color,
  node_color_light,
  save_file = F,
  filename = "outfile.rmd",
  cores = 6
){
  data_mat <- t(assay(SE))
  sample_data <- data.frame(colData(SE))
  mol_data <- data.frame(rowData(SE))
  R <- initialize(data_mat, sample_data, mol_data, cores) %>% 
    find_sig_clusts(phenotype, confounders, cores, node_color, node_color_light, nrand = 1000)
  R$phenotypes = phenotype
  if(save_file){
    save(R, file = filename)
  }
  R
}


##################################################################
