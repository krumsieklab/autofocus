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
#' @param use_wgcna Whether or not to use WGNA for hiearchical clustering. Default: F.
#' @param adj_method Which p-value adjustment method to use. One of c("bonferroni", "fdr").
#'
#' @return R struct with the data matrix, correlation matrix,
#' sample and molecular annotations, hierarchical structure, cluster membership
#' node order within the hierarchical structure and a color vector for the platforms
#'
#' @author AS
#'
#' @importFrom dendextend nnodes
#' @import parallel
#' @import foreach
#' @importFrom doParallel registerDoParallel
#'
#' @export
initialize_R <- function(data.matrix,
                         sample.data,
                         mol.data,
                         confounders,
                         phenotype,
                         cores = 4,
                         use_wgcna = F,
                         adj_method = c("bonferroni", "fdr")
){

  adj_method = match.arg(adj_method)

  # check that data.matrix is input the correct way (samples in rows and features in columns)
  # crash if not
  df_rows <- nrow(data.matrix)
  df_cols <- ncol(data.matrix)
  nsamp <- nrow(sample.data)
  nmol <- nrow(mol.data)
  # check that number of rows in data.matrix is equal to number of rows in sample.data
  if(df_rows != nsamp){
    if(df_rows == nmol){
      stop("Number of rows of sample.data do not match number of data.matrix rows. Did you forget
      to transpose SE assay?")
    }else{
      stop("Number of rows of data.matrix (samples) do not match rows in sample.matrix.")
    }
  }
  # check that number of columns in data.matrix is equal to number of rows in mol.data
  if(df_cols != nmol){
    if(df_cols == nsamp){
      stop("Number of rows of mol.data do not match number of data.matrix columns. Did you forget
      to transpose SE assay?")
    }else{
      stop("Number of columns of data.matrix (features) do not match rows in mol.matrix.")
    }
  }

  # Organize input data
  R <- list(data = scale(data.matrix), samples = sample.data, annos = mol.data)

  # run WGCNA?
  if(use_wgcna){
    dissTom <- get_wgcna_dist_metric(data.matrix, cores)
    R$HCL <- hclust(as.dist(dissTom),method="average")
  }else{
    # Calculate correlations between biomolecules
    R$C <- stats::cor(data.matrix, use="pairwise.complete.obs")

    # Calculate hierarchical structure, order of nodes
    R$HCL <- R %>% get_dendro(method = "average", cores)
  }

  R$dend_data <- try(ggdendro::dendro_data(as.dendrogram(R$HCL), type = "rectangle"), silent = T)
  if(class(R$dend_data)=="try-error"){
    if(R$dend_data == "Error : node stack overflow\n"){
      warning("Hierarchical tree exceeds maximum recursion depth and cannot be drawn with the
              ggdendro::dendro_data function. This will affect visualization of the dendrogram
              in the autofocus shiny app.")
    }else{
      stop(glue::glue("The following error was encountered running ggdendro::dendro_data:\n", R$dend_data))
    }
  }

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
    R$clust_info$densities <- get_sig_child_density(R, phenotype, confounders, method=adj_method)

    R$clust_info$mgm_capable <- (num_samples >= R$clust_info$Size) & (R$clust_info$densities >= 0.5) & (R$clust_info$Size>1)

    if(cores > 1){
      registerDoParallel(cores=cores)
      R$graphs <- foreach(i=1:nnodes(R$HCL)) %dopar% {get_edges_linear(i, R, phenotype, confounders)}
    }else{
      R$graphs <- lapply(1:nnodes(R$HCL), function(i) get_edges_linear(i,R, phenotype, confounders))
    }

  }else if(length(phenotype)==2){
    R$clust_info$pheno1_densities <- get_sig_child_density(R, phenotype[1], confounders, method=adj_method)
    R$clust_info$pheno1_mgm_capable <- (num_samples >= R$clust_info$Size) & (R$clust_info$pheno1_densities >= 0.5) & (R$clust_info$Size>1)

    R$clust_info$pheno2_densities <- get_sig_child_density(R, phenotype[2], confounders, method=adj_method)
    R$clust_info$pheno2_mgm_capable <- (num_samples >= R$clust_info$Size) & (R$clust_info$pheno2_densities >= 0.5) & (R$clust_info$Size>1)

    if(cores > 1){
      registerDoParallel(cores=cores)
      R$graphs <- foreach(i=1:nnodes(R$HCL)) %dopar% {get_edges_linear(i, R, phenotype, confounders)}
    }else{
      R$graphs <- lapply(1:nnodes(R$HCL), function(i) get_edges_linear(i,R, phenotype, confounders))
    }


  }else{
    stop("autofocus only supports viewing up to two phenotypes.")
  }
  to_remove <- c("data", "dist","C","order")
  R[!(names(R) %in% to_remove)]

}

