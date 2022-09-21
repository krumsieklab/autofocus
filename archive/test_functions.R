bind_SE_with_NA <- function(platforms, sample_id){
  
  # Set the column names of each platform to the universal sample identifier
  platforms <- lapply(platforms, function(i) {
    colnames(i)<-colData(i)[[sample_id]]
    i}
  )
  
  # Bind the assay data, order
  full_assay <- lapply(platforms, function(i) data.frame(assay(i), check.names=FALSE)) %>% bind_rows()
  
  columns_df <- lapply(platforms, function(i) data.frame(colData(i)))
  
  full_colData <- Reduce(function(x, y) merge(x, y, by=colnames(columns_df[[1]]), all=T), columns_df)
  rownames(full_colData) <- full_colData$SampleId
  full_colData <- full_colData[colnames(full_assay),]
  
  full_rowData <- matchRowClasses(platforms)
  
  SummarizedExperiment(assays = full_assay, colData = full_colData, rowData = full_rowData)
}



matchRowClasses <-  function(platforms){
  
  rows_df <- lapply(platforms, function(i) data.frame(rowData(i)))
  
  row_names_types <- lapply(rows_df, function(x) {
    data.frame(Name = names(x), Type = sapply(x, typeof))
  }) %>% 
    bind_rows()  %>%
    group_by(Name) %>%
    slice(1) %>%
    ungroup() %>% 
    droplevels()
  
  
  for(x in 1:length(rows_df)){
    for (i in 1:ncol(rows_df[[x]])){
      class(rows_df[[x]][,i])<- as.character(row_names_types$Type[which(row_names_types$Name == colnames(rows_df[[x]])[i])])
    }
    
  } 
  
  rows_df %>% bind_rows()
  
}


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
                             Size = c(mapply(function(i)length(i),R$clusts),rep(1, nleaves(R$HCL))))
  
  R$clust_info$densities <- get_sig_child_density(R, phenotype, confounders)
  
  R$clust_info$mgm_capable <- (num_samples >= R$clust_info$Size) & (R$clust_info$densities >= 0.5) & (R$clust_info$Size>1)
  R$graphs <- lapply(1:nnodes(R$HCL), function(i) get_edges_linear(i,R, phenotype, confounders))
  to_remove <- c("data", "dist","C","order")
  R[!(names(R) %in% to_remove)] 
}
