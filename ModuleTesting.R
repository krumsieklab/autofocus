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
    
    R$BIC <- foreach(i = 1:nnodes(R$HCL)) %dopar% {scoring_func_wrapper(i, 
                                                                        R, 
                                                                        as.matrix(R$samples[[phenotype]]), 
                                                                        confounders, 
                                                                        score_method,
                                                                        return_BIC = T) }
    
    R$pvals <- p_adjust_wrapper(unlist(allpvals),
                                inds,
                                R,
                                as.matrix(R$samples[[phenotype]]),
                                confounders,
                                score_method = score_method,
                                adjust_method)
  
    # determine significant nodes to be colored
    signif <- which(R$pvals<0.05)
    
    R$colors <- mapply(function(i) get_node_color(R, i, signif, node_color, node_color_light), 1:nnodes(R$HCL))
  
    to_remove <- c("data", "dist")
    R[!(names(R) %in% to_remove)]  
}
