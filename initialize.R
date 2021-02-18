### internal functions

#source("~/Documents/JanLab/git/autofocus/Backend_Modules.R")

#### Parse Input ####
#' input
#' 
#' Parses input data
#'
#' @param mol_sheets List of excel file names containing molecular data
#' @param sample_sheet Excel file name containing sample data
#' @param id_col string, name of column in sample_sheet that represents identifiers of samples, used in all mol_sheets as well
#' @param phenotype string, name of column containing phenotype to be tested
#' @param confounders vector of strings, names of columns in sample_sheet that contain confounders
#' @param cores Number of cores to be used in algorithm
#' 
#' @return R struct with the data matrix, correlation matrix, 
#' sample and molecular annotations, hierarchical structure, cluster membership
#' node order within the hierarchical structure and a color vector for the platforms
input <- function(
  mol_sheets,
  sample_sheet,
  id_col,
  phenotype,
  confounders,
  cores
){
  # TODO: read in sheets. Each mol_sheet should have
  # two excel sheets, one is data (preprocessed), one is molecular annotations
  # TODO: checks
  #         no NAs
  #         IDs match
  #         Dimensions match
  #         as.numeric
  # Add platform column to molecular data
  
  
}

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

initialize_R <- function(
  data.matrix, # Matrix of raw data (samples * nodes)
  sample.data, # Matrix of sample annotations
  mol.data, # Matrix of node annotations
  cores
){
  
  # Organize input data
  R <- list(data = data.matrix, samples = sample.data, annos = mol.data)
  # Calculate correlations between biomolecules
  R$C <- cor(data.matrix)
  
  # Calculate hierarchical structure, order of nodes
  R$HCL <- R %>% get_dendro(method = "average", cores)
  R$dend_data <- ggdendro::dendro_data(as.dendrogram(R$HCL), type = "rectangle")
  R$order <- get_dend_indices(R, dim(R$HCL$merge)[1], c())
  
  # Get all clusters in the hierarchical tree
  R$clusts <- lapply(1:dim(R$HCL$merge)[1], function(x) get_members(R, x))
  dend_xy <- R$HCL %>% as.dendrogram %>%get_nodes_xy()
  
  # Get coordinate information
  parents <- unlist(lapply(1:nnodes(R$HCL), function(i) get_parent(R,i)))
  coord_x <- unlist(lapply(1:nnodes(R$HCL), function(i) round(dend_xy[which(R$order==i),1], digits = 10)))
  coord_y <- unlist(lapply(1:nnodes(R$HCL), function(i) round(dend_xy[which(R$order==i),2],digits= 10)))
  
  R$clust_info <- data.frame(ClusterID=1:nnodes(R$HCL), 
                             Parent = parents, 
                             Coord_X=coord_x, 
                             Coord_Y=coord_y,
                             Size = c(mapply(function(i)length(i),R$clusts),rep(1, nleaves(R$HCL))))
  
  R
}

#### Get parent of a cluster ####
#' get_parent
#' 
#' Recursively collects the indices of 
#' leaves that are descendents of a module
#'
#' @param R R struct
#' @param i index of module we are finding the members of
#' @param cl vector to collect node indices
#' 
#' @return vector of indices of leaves in a module 

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

