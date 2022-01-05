### internal functions
source("Backend_Modules.R")

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
#' @param confounders Vector of names of potential confounders (names must be in sample.data)
#' @param phenotype Name of phenotype of interest, must be in sample.data
#' @param cores Number of cores to run on
#' 
#' @return R struct with the data matrix, correlation matrix, 
#' sample and molecular annotations, hierarchical structure, cluster membership
#' node order within the hierarchical structure and a color vector for the platforms

initialize_R <- function(
  data.matrix, # Matrix of raw data (samples * nodes)
  sample.data, # Matrix of sample annotations
  mol.data, # Matrix of node annotations
  confounders, # Vector of confounder names (must be found in sample data)
  phenotype, # Name of column containing disease phenotype of interest
  cores = 4
){
  
  # Organize input data
  R <- list(data = scale(data.matrix), samples = sample.data, annos = mol.data)
  # Calculate correlations between biomolecules
  R$C <- stats::cor(data.matrix, use="pairwise.complete.obs")
  
  # Calculate hierarchical structure, order of nodes
  R$HCL <- R %>% get_dendro(method = "average", cores)
  R$dend_data <- ggdendro::dendro_data(as.dendrogram(R$HCL), type = "rectangle")
  R$order <- get_dend_indices(R, dim(R$HCL$merge)[1], c())
  
  # Get all clusters in the hierarchical tree
  R$clusts <- lapply(1:dim(R$HCL$merge)[1], function(i) {get_members(R, i)})
  dend_xy <- R$HCL %>% as.dendrogram %>%get_nodes_xy()
  
  # Get coordinate information
  parents <- lapply(1:nnodes(R$HCL), function(i) {get_parent(R,i)}) %>% unlist()
  coord_x <-lapply(1:nnodes(R$HCL), function(i) {round(dend_xy[which(R$order==i),1], digits = 10)}) %>% unlist()
  coord_y <- lapply(1:nnodes(R$HCL), function(i){round(dend_xy[which(R$order==i),2],digits= 10)}) %>% unlist()
  
  num_samples <- lapply(1:nnodes(R$HCL), function (i) {
    if (i > dim(R$HCL$merge)[1]){
      return(sum(!is.na(R$data[,(i-dim(R$HCL$merge)[1])])))
    }
    else{
      members <-  R$clusts[i][[1]]
      data = R$data[,members]
      return(1*(apply(is.na(data),1,sum)==0) %>% sum())
    }
  }) %>% unlist()

  R$clust_info <- data.frame(ClusterID=1:nnodes(R$HCL),
                             Parent = parents,
                             Coord_X=coord_x,
                             Coord_Y=coord_y,
                             Size = c(mapply(function(i)length(i),R$clusts),rep(1, nleaves(R$HCL))))
  
  R$clust_info$densities <- get_sig_child_density(R, phenotype, confounders)
  
  R$clust_info$mgm_capable <- (num_samples >= R$clust_info$Size) & (R$clust_info$densities > 0.5) & (R$clust_info$Size>1)
  R$graphs <- lapply(1:nnodes(R$HCL), function(i) get_edges_linear(i,R, phenotype, confounders))
  to_remove <- c("data", "dist","C","order")
  R[!(names(R) %in% to_remove)] 
}

