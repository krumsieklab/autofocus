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
  cores = 1, 
  node_color,
  node_color_light,
  nrand = 1000)
  {
    if (cores>1) doParallel::registerDoParallel(cores=cores)
  
    #### find all valid clusters and score them ----
    PHEN <- R$samples[[phenotype]]
    ## calculate p-values
    all_nodes <- nnodes(R$HCL)
    
    print("Calculating pvals")
    
    allpvals <- foreach(r = 1:all_nodes) %dopar% {module.pval(r, R, PHEN, confounders, max_size = 50) }
    inds <- which(!is.na(allpvals))
    
    print("Calculating BIC")
    
    R$BIC <- foreach(r = 1:all_nodes) %dopar% {module.pval(r, R, PHEN, confounders, TRUE, max_size = 50) }
    
    
    #### WY p-values ----
    ## -> randomize outcome, run all tests, record smallest p-value of each iterations
    print("Calculating WY pvals")
    
    wy.null <- foreach(r = 1:nrand) %dopar% {  min(sapply(inds, module.pval, R, sample(PHEN), confounders, max_size = 50)) } %>% unlist()
    
    # NAs can occur (very rarely), cut them out
    cut <- which(is.na(wy.null))
    if (length(cut)>0) {
      print(sprintf("NAs found in random samples: %s", paste0(cut, collapse = ',')))
      wy.null <- wy.null[-cut]
      nrand <- length(wy.null)
    }
    
    ## empirical WY adjusted p-value generation
    # -> count how many of the WY null values are actually lower than the observed p-value
    R$pvals <- sapply(allpvals, function(pval){sum(wy.null<pval)}) / nrand
  
    # determine significant nodes to be colored
    signif <- which(R$pvals<0.05)
    
    R$colors <- mapply(function(i) get_node_color(R, i, signif, node_color, node_color_light), 1:all_nodes)
    R$labels <- mapply(function(i) get_node_label(R, i), 1:all_nodes)
  
    to_remove <- c("data", "dist")
    R[!(names(R) %in% to_remove)]  
}


#' module.pval
#' 
#' Takes as input the hierarchical structure, a module within that structure and molecular and phenotype information
#' and finds the p-value or BIC of the module
#'
#' @param i Index of parent of the module
#' @param R R struct
#' @param phenotype_vec Vector of phenotype values 
#' @param confounders Vector of confounder names as strings
#' @param return_BIC Return BIC instead of p-value?
#' @param max_size Maximum size a module can be
#' 
#' @return p-value or BIC of module


module.pval <- function(
  i,
  R,
  phenotype_vec, 
  confounders, 
  return_BIC = FALSE, 
  max_size
){
  X <- R$data
  full_data <- cbind(X, phenotype_vec, as.matrix(R$samples[,confounders]))
  me <- R$HCL$merge
  if (i <= dim(me)[1]){
    members <- R$clusts[i][[1]]
    if (length(members)>max_size){
      return(NA)
    }
  }else{
    members <- (i - dim(me)[1])
  }
  
  if (length(confounders)==0){
    # test model
    model.test <- lm(phenotype_vec ~ full_data[,members] + 1) 
    # base model
    model.base <- lm(phenotype_vec ~ 1)
  }else{
    # test model
    model.test <- lm(phenotype_vec ~ full_data[,members] + full_data[,confounders]+1) 
    # base model
    model.base <- lm(phenotype_vec ~ full_data[,confounders]+1) 
  }
  # BIC calculation
  if (return_BIC){
    return(BIC(model.test))
  }
  # p-val of model differences
  anova(model.test, model.base)$`Pr(>F)`[2]
}
