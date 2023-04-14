#' Initialize autofocus data struct
#'
#' Takes as input a data matrix, annotations on the samples and annotations on the molecules
#' and creates a list structure, R with the hierarchical clustering structure and cluster informations.
#'
#' @param data.matrix Matrix with samples as rows and biomolecules as columns.
#' @param sample.data Matrix of annotation data of samples.
#' @param mol.data Matrix of annotation data of biomolecules.
#' @param confounders Vector of names of potential confounders (names must be in sample.data).
#' @param phenotype Name of column containing phenotype of interest, must be in sample.data.
#' @param cores Number of cores to run on. Default: 4.
#'
#' @return R struct with the data matrix, correlation matrix,
#' sample and molecular annotations, hierarchical structure, cluster membership
#' node order within the hierarchical structure and a color vector for the platforms
#'
#' @author AS
#'
#' @importFrom dendextend nnodes
#'
#' @export
initialize_R <- function(data.matrix,
                         sample.data,
                         mol.data,
                         confounders,
                         phenotype,
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
  dend_xy <- R$HCL %>% as.dendrogram %>% dendextend::get_nodes_xy()

  # Get coordinate information
  parents <- lapply(1:nnodes(R$HCL), function(i) {get_parent(R,i)}) %>% unlist()
  coord_x <-lapply(1:nnodes(R$HCL), function(i) {round(dend_xy[which(R$order==i),1], digits = 10)}) %>% unlist()
  coord_y <- lapply(1:nnodes(R$HCL), function(i){round(dend_xy[which(R$order==i),2],digits= 10)}) %>% unlist()

  # Get the number of full samples in each cluster
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

  # Build data.frame containing data for each cluster
  R$clust_info <- data.frame(ClusterID=1:nnodes(R$HCL),
                             Parent = parents,
                             Coord_X=coord_x,
                             Coord_Y=coord_y,
                             Size = c(mapply(function(i)length(i),R$clusts),rep(1, dendextend::nleaves(R$HCL))))

  if(length(phenotype)==1){
    R$clust_info$densities <- get_sig_child_density(R, phenotype, confounders)

    R$clust_info$mgm_capable <- (num_samples >= R$clust_info$Size) & (R$clust_info$densities >= 0.5) & (R$clust_info$Size>1)
    R$graphs <- lapply(1:nnodes(R$HCL), function(i) get_edges_linear(i,R, phenotype, confounders))

  }else if(length(phenotype)==2){
    R$clust_info$pheno1_densities <- get_sig_child_density(R, phenotype[1], confounders)
    R$clust_info$pheno1_mgm_capable <- (num_samples >= R$clust_info$Size) & (R$clust_info$pheno1_densities >= 0.5) & (R$clust_info$Size>1)

    R$clust_info$pheno2_densities <- get_sig_child_density(R, phenotype[2], confounders)
    R$clust_info$pheno2_mgm_capable <- (num_samples >= R$clust_info$Size) & (R$clust_info$pheno2_densities >= 0.5) & (R$clust_info$Size>1)

    R$graphs <- lapply(1:nnodes(R$HCL), function(i) get_edges_linear(i,R, phenotype, confounders))

  }else{
    stop("autofocus only supports viewing up to two phenotypes.")
  }
  to_remove <- c("data", "dist","C","order")
  R[!(names(R) %in% to_remove)]

}

