#### File that takes as input the R structure and performs disease association testing ###
source(codes.makepath("autofocus/internal_funcs.R"))

### Main function, returns R structure with nodes labelled by significance----------------

find_sig_clusts <- function(R, outcome, confounders, cores = 1, nrand = 1000){
  hc <- R$HCL
  if (cores>1) doParallel::registerDoParallel(cores=cores)
  PHEN <- R$samples[[outcome]]
  X <- R$data
  full_data <- cbind(X, PHEN, as.matrix(R$samples[,confounders]))
  
  #### find all valid clusters and score them ----

  ## calculate p-values
  all_nodes <- nnodes(hc)
  
  print("Calculating pvals")
  
  allpvals <- foreach(r = 1:all_nodes) %dopar% {module.pval(r, full_data, R, PHEN, confounders, max_size = 50) }
  inds <- which(!is.na(allpvals))
  
  print("Calculating BIC")
  
  R$BIC <- foreach(r = 1:all_nodes) %dopar% {module.pval(r, full_data, R, PHEN, confounders, TRUE, max_size = 50) }
  
  
  #### WY p-values ----
  ## -> randomize outcome, run all tests, record smallest p-value of each iterations
  print("Calculating WY pvals")
  
  wy.null <- foreach(r = 1:nrand) %dopar% {  min(sapply(inds, module.pval, full_data, R, sample(PHEN), confounders, max_size = 50)) } %>% unlist()
  
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
  
  R$colors <- mapply(function(i) get_node_color(R, i, signif), 1:all_nodes)
  R$labels <- mapply(function(i) get_node_label(R, i), 1:all_nodes)

  to_remove <- c("data", "dist")
  R[!(names(R) %in% to_remove)]  
}


### Calculates the p-value or BIC of a module ----------------
module.pval <- function(i, full_data, R, phenotype, confounders, return_BIC = FALSE, max_size){
  me <- R$HCL$merge
  if (i <= dim(me)[1]){
    members <- R$clusts[i][[1]]
    if (length(members)>max_size){
      return(NA)
    }
  }
  else{
    members <- (i - dim(me)[1])
  }
  
  if (length(confounders)==0){
    # test model
    model.test <- lm(phenotype ~ full_data[,members] + 1) 
    # base model
    model.base <- lm(phenotype ~ 1)
  }
  else{
    # test model
    model.test <- lm(phenotype ~ full_data[,members] + full_data[,confounders]+1) 
    # base model
    model.base <- lm(phenotype ~ full_data[,confounders]+1) 
  }
  # BIC calculation
  if (return_BIC){
    return(BIC(model.test))
  }
  # p-val of model differences
  anova(model.test, model.base)$`Pr(>F)`[2]
}
