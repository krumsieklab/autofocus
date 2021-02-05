#### File that takes as input the R structure and performs disease association testing ###

#' find_sig_clusts
#' 
#' Takes as input the hierarchical structure and molecular and phenotype information
#' and performs significance testing over the modules of the tree
#'
#' @param R R struct 
#' @param phenotype Phenotype name 
#' @param confounders Vector of confounder names as strings
#' @param cores Number of cores to use for parallelization
#' @param nrand Number of randomizations for WY null distribution
#' 
#' @return R struct with nodes labeled and colored based on BIC and p-value
#'

find_sig_clusts <- function(
  R, 
  phenotype, 
  confounders,
  score_method = "lm",
  adjust_method = "wy",
  cores = 1, 
  node_color = "green",
  node_color_light = "lightgreen")
  {
    if (cores>1) doParallel::registerDoParallel(cores=cores)
   
   ## calculate p-values
    
    print("Calculating pvals")
    
    allpvals <- foreach(i = 1:nnodes(R$HCL)) %dopar% {scoring_func_wrapper(i, 
                                                                           R, 
                                                                           R$samples[[phenotype]], 
                                                                           confounders, 
                                                                           score_method) }
    inds <- which(!is.na(allpvals))

    print("Calculating BIC")
    
    R$clust_info$BIC <- foreach(i = 1:nnodes(R$HCL)) %dopar% {scoring_func_wrapper(i, 
                                                                        R, 
                                                                        as.matrix(R$samples[[phenotype]]), 
                                                                        confounders, 
                                                                        score_method,
                                                                        return_BIC = T) } %>% unlist() %>% round(digits=5)
    
    R$clust_info$pvals <- p_adjust_wrapper(unlist(allpvals),
                                inds,
                                R,
                                as.matrix(R$samples[[phenotype]]),
                                confounders,
                                score_method = score_method,
                                adjust_method) %>% round(digits=5)
  
    # determine significant nodes to be colored
    signif <- which(R$clust_info$pvals<0.05)
    print(signif)
    R$clust_info$colors <- mapply(function(i) get_node_color(R, i, signif, node_color, node_color_light), 1:nnodes(R$HCL))
  
    to_remove <- c("data", "dist","C")
    R[!(names(R) %in% to_remove)]  
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
      if (R$clust_info$pvals[child1] > 0.05 & R$clust_info$pvals[child2] > 0.05) return(node_color)
      
      ## Children are significant
      if(R$clust_info$pvals[child1] < 0.05 & R$clust_info$pvals[child2] < 0.05){
        
        if (R$clust_info$BIC[[i]] < R$clust_info$BIC[[child1]] & R$clust_info$BIC[[i]]<R$clust_info$BIC[[child2]]) return(node_color) else return(node_color_light)
      }
      
      ## One child is significant
      else {
        sig_child <- if (R$clust_info$pvals[child1] < R$clust_info$pvals[child2]) child1 else child2
        
        if(R$clust_info$BIC[[i]] < R$clust_info$BIC[[sig_child]]) return(node_color) else return(node_color_light)
      }
    }
  }
  else{
    return("black")
  }
}
