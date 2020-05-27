#### init ----
source(codes.makepath("autofocus/internal_funcs.R"))


get_node_color <- function(R, i, signif){
  hc <- R$HCL
  internal_nodes <- dim(hc$merge)[1]
  if (i %in% signif[!is.na(signif)]){
    if (i > (internal_nodes)) return('blue')
    else {
      
      children <- hc$merge[i,]
      child1 <- if (children[1]<0) (internal_nodes + abs(children[1])) else children[1]
      child2 <- if (children[2]<0) (internal_nodes + abs(children[2])) else children[2]
      
      ## Children aren't significant, combination is
      if (R$pvals[child1] > 0.05 & R$pvals[child2] > 0.05) return("purple")
      
      ## Children are significant
      if(R$pvals[child1] < 0.05 & R$pvals[child2] < 0.05){
        
        if (R$BIC[i] < R$BIC[child1] & R$BIC[i]<R$BIC[child2]) return("purple") else return("orange")
      }
      
      ## One child is significant
      else {
        sig_child <- if (R$pvals[child1] < R$pvals[child2]) child1 else child2
        
        if(R$BIC[i] < R$BIC[sig_child]) return("purple") else return("orange")
      }
    }
  }
  else{
    return("black")
  }
}

get_node_label <- function(R, i){
  internal_nodes <- dim(R$HCL$merge)[1]
  if (i > (internal_nodes)) return(R$HCL$labels[(i - internal_nodes)])
  else return(paste("BIC: ", round( R$BIC[i], digits = 3), ", p-value: ", R$pvals[i]))
}


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
  # p-val of model differences
  if (return_BIC){
    return(BIC(model.test))
  }
  anova(model.test, model.base)$`Pr(>F)`[2]
}

find_sig_clusts <- function(R, outcome, confounders, cores = 1){
  hc <- R$HCL
  
  PHEN <- R$samples[[outcome]]
  X <- R$data
  #### scoring function definition ----
  full_data <- cbind(X, PHEN, as.matrix(R$samples[,confounders]))
  
  #### find all valid clusters and score them ----

  ## calculate p-values
  internal_nodes <- dim(hc$merge)[1]
  leaves <- nleaves(hc)
  all_nodes <- nnodes(hc)
  print("Calculating pvals")
  allpvals <- sapply(1:all_nodes, module.pval, full_data, R, PHEN, confounders, max_size = 50); 
  print("Calculating BIC")
  R$BIC <- sapply(1:all_nodes, module.pval, full_data, R, PHEN, confounders, max_size = 50, TRUE);
  
  inds <- which(!is.na(allpvals))
  
  #### WY p-values ----
  ## -> randomize outcome, run all tests, record smallest p-value of each iterations
  
  # server timing, with old version, not all nodes (singltons missing, real time will be about double)
  
  # - nrand=100000, 40 cores, 135184.517 sec elapsed (37.5h)
  
  # parameters
  nrand <- 50
  mc.cores <- cores
  # safe parallelization
  
  # run
  
  if (mc.cores>1) doParallel::registerDoParallel(cores=mc.cores)
  wy.null <- foreach(r = 1:nrand) %dopar% {  min(sapply(inds, module.pval, full_data, R, sample(PHEN), confounders, max_size = 50)) } %>% unlist()
  
  # NAs can occur (very rarely), cut them out
  cut <- which(is.na(wy.null))
  if (length(cut)>0) {
    print(sprintf("NAs found in random samples: %s", paste0(cut, collapse = ',')))
    wy.null <- wy.null[-cut]
    nrand <- length(wy.null)
    # plot null distribution
    hist(log10(wy.null))
    plot(log10(wy.null))
  }
  
  ## empirical WY adjusted p-value generation
  # -> count how many of the WY null values are actually lower than the observed p-value
  R$pvals <- sapply(allpvals, function(pval){sum(wy.null<pval)}) / nrand

  # determine significant nodes to be colored
  print("Getting node colors")
  signif <- which(R$pvals<0.05)
  
  R$colors <- mapply(function(i) get_node_color(R, i, signif), 1:all_nodes)
  R$labels <- mapply(function(i) get_node_label(R, i), 1:all_nodes)
  
  R<- R %>% make_bar_plots()
  to_remove <- c("data", "dist")
  R[!(names(R) %in% to_remove)]  
}