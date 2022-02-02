## Backend Modules ###

suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(mgm))

library(RhpcBLASctl)
blas_set_num_threads(1)


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


## Get dendrogram ##
#' get_dendro
#' 
#' Hierarchically clusters data in R
#'
#' @param R R struct
#' @param method Hierarchical clustering method. If NA (default), get optimal choice of all clustering methods
#' @param cores Number of cores to be used in this step
#' 
#' @return hclust object of dendrogram
#' 
get_dendro <- function(
  R,
  method = NA,
  cores = cores
){
  
  # Calculate distance metric
  dist_mat <- abs_cor_dist(R)
  
  # Hclust using defined method if not NA
  if (!(is.na(method))){
    return(hclust(dist_mat, method = method))
  }
  
  # Hclust using search over all clustering methods
  else {
    return(get_cluster_method(dist_mat, cores = cores))
  }
}


## Distance Metric ##
#' abs_cor_dist
#' 
#' Calculates distance matrix of data based on correlation value
#'
#' @param R R struct
#' 
#' @return Distance matrix made up of 1- the absolution value of the correlation matrix
#' 
abs_cor_dist <- function(
  R
){
  as.dist((1-abs(R$C))) # abs correlation
}


## Get dendrogram from best clustering method ##
#' get_cluster_method
#' 
#' Searches through all hierarchical clustering methods for one that
#' results in a dendrogram with highest cophenetic correlation to original data matrix
#'
#' @param dist_mat distance matrix
#' @param cores Number of cores to be used in this step
#' 
#' @return hclust object of dendrogram
#' 
get_cluster_method <- function(
  dist_mat, 
  cores = 1
){
  
  # Test values for agnes clustering function
  test_vals <- seq(0,1,by=0.05)
  
  #For each test value, return correlation coefficient between distance matrix
  # and cophenetic values of the resulting dendrogram using the flexible parameter
  flex_known_cors <- unlist(mclapply(test_vals, function(i){
    hc = as.hclust(cluster::agnes(dist_mat, diss = T,method = "flexible", par.method = i))
    cor(as.vector(as.dist(cophenetic(hc))), as.vector(dist_mat), method = "spearman")
  }, mc.cores = cores))
  
  flex_known_data <- data.frame(flex_known_cors, test_vals,clustering_type = rep("flexible", length(flex_known_cors)))
  
  #For each test value, return correlation coefficient between distance matrix
  # and cophenetic values of the resulting dendrogram using the gaverage parameter
  gav_known_cors <- unlist(mclapply(test_vals, function(i){
    hc = as.hclust(cluster::agnes(dist_mat, diss = T,method = "gaverage", par.method = i))
    cor(as.vector(as.dist(cophenetic(hc))), as.vector(dist_mat), method = "spearman")
  }, mc.cores = cores))
  gav_known_data <- data.frame(gav_known_cors, test_vals, clustering_type = rep("gaverage", length(gav_known_cors)))
  
  # Plot the values for each method
  known_data <- rbind(flex_known_data, gav_known_data)
  p<- ggplot(known_data) +geom_point(aes(x=test_vals, y=known_cors))
  print(p + 
          facet_grid(cols = vars(clustering_type)) +
          ggtitle("Correlation between linkage method and distance matrix") +
          xlab("value of alpha") + ylab("Cophenetic correlation") + 
          theme_minimal() + 
          theme(strip.background = element_rect(fill = "wheat", color = "black"),
                panel.background = element_rect(colour = "black") ) + 
          ggsci::scale_color_aaas())
  
  # Return the highest correlated dendrogram 
  if (max(flex_known_cors)>max(gav_known_cors)){
    best_flex_ind <- which(flex_known_cors == max(flex_known_cors))
    return(as.hclust(cluster::agnes(dist_mat, 
                                    diss = T,
                                    method = "flexible", 
                                    par.method = test_vals[best_flex_ind])))
  }
  else{
    best_gav_ind <- which(gav_known_cors == max(gav_known_cors))
    return(as.hclust(cluster::agnes(dist_mat, 
                                    diss = T,
                                    method = "gaverage", 
                                    par.method = test_vals[best_gav_ind])))
  }
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

#### Get proportion/density of significant children ####
#' get_sig_child_density
#' 
#' Finds the proportion of sigificant leaves at each
#' parent node of the dendrogram 
#'
#' @param R R struct
#' @param phenotype disease pheotype 
#' @param confounders vectors of confounder names 
#' 
#' @return Returns the densities and number of significant descendants of each node
#' 
get_sig_child_density <- function(R, phenotype, confounders){
  sig_kids<-lapply(1:nleaves(R$HCL), function(i) {
    data <- data.frame(pheno=R$samples[[phenotype]], molecule = R$data[,i], R$samples[,confounders])
    reg <- lm(pheno ~ ., data = data)
    summary(reg)$coefficients["molecule", "Pr(>|t|)"]
    
  }) %>% unlist %>%  p.adjust(method="bonferroni")
  
  sig_kids<-1*(sig_kids<0.05)
  
  densities <- lapply(1:dim(R$HCL$merge)[1], function(i) {
    members <-  R$clusts[i][[1]]
    sum(sig_kids[members])/length(members)
  }) %>% unlist()
  c(densities, sig_kids)
}

#### Get direct edges between disease phenotype and module members ####
#' get_edges_linear
#' 
#' Uses mixed graphical models to create a graph
#' of direct edges between module members, confounders,
#' and disease phenotype
#'
#' @param i index of parent node of cluster
#' @param R R struct
#' @param phenotype disease pheotype 
#' @param confounders vectors of confounder names 
#' 
#' @return Returns a graph object
#'
get_edges_linear <- function(i, R, phenotype, confounders){

  if (!(R$clust_info$mgm_capable[i])){
    return (NA)
  } 
  
  else{
    members <-  R$clusts[i][[1]]
    
    # Get variable type of phenotype and confounders
    # Either categorical ("c") or gaussian ("g")
    # If categorical, get number of categories
    
    #phenotype_type <- get_variable_type(R$samples[[phenotype]])
    #phenotype_cat <- get_num_categories(R$samples[[phenotype]])
    
    confounder_data <- R$samples[,confounders]
    #confounder_types <- apply(confounder_data, 2, get_variable_type) %>% as.vector()
    #confounder_cat <- apply(confounder_data, 2, get_num_categories) %>% as.vector()
    
    #Combine all data with confounders
    full_data<-data.frame(R$data[,members], confounder_data)
    names <- c(phenotype, names(full_data))
    full_data <-as.matrix(as.data.frame(lapply(full_data, as.numeric)))
    # Remove samples with NAs
    no_nas = apply(is.na(full_data),1,sum) == 0
    data_no_na <- cbind(R$samples[[phenotype]][no_nas], full_data[no_nas,])
    colnames(data_no_na)<- names
    
    types <- apply(data_no_na, 2, get_variable_type) %>% as.vector()
    levels <- apply(data_no_na, 2, get_num_categories) %>% as.vector()
    
    if('c' %in% types) {
      ind_cat <- which(types == 'c')
      to_remove<-c()
      for(i in ind_cat) {
        l_frqu <- table(data_no_na[,i])
        cats_to_remove <- names(l_frqu)[which(l_frqu<=1)]
        if(length(cats_to_remove)>0){
          levels[i]<-levels[i] - length(cats_to_remove)
          to_remove<-c(to_remove, which(data_no_na[,i]%in% cats_to_remove))
          }
        } # this does not catch the case where one category is not present at all; but this is catched by comparing specified levels and real levels
    }
    if(length(to_remove)>0){
      data_no_na <- data_no_na[-to_remove,]
    }
    # Find direct connections with phenotype in mgm
    fit_mgm <- mgm(data_no_na, type=types, level=levels)
    adj_mat <- 1* (fit_mgm$pairwise$wadj!=0)
    
    # Create graph
    colnames(adj_mat)<- c(phenotype, colnames(R$data)[members], confounders)
    G<-graph_from_adjacency_matrix(adj_mat, add.rownames = T)
    V(G)$NodeType <- c("phenotype", rep("analyte",length(members)), rep("confounder", length(confounders)))
    V(G)$Driver <- V(G)$name %in% neighbors(G,phenotype)$name[neighbors(G,phenotype)$NodeType =="analyte"]
    V(G)$color <- c("blue", rep("red",length(members)), rep("grey", length(confounders)))
    G
  }
}

#### Determine whether a variable is discrete or continuous ####
#' get_variable_type
#' 
#' Determines whether a variable is discrete or continuous based
#' on the ratio of the number of unique entries to number of entries
#'
#' @param variable_vec The vector to determine type
#' @param cutoff The cutoff of the unique entries to length ratio 
#' 
#' @return Returns "c" if discrete categorial or "g" if continuous or gaussian
#'
get_variable_type <- function(variable_vec, cutoff= 0.05){
  ifelse((length(unique(variable_vec))/length(variable_vec)) <= cutoff, "c", "g")
}

#### Determine number of categories in a variable ####
#' get_num_categories
#' 
#' If a variable is determinned to be categorical, the number of categories
#' Else returns 1
#'
#' @param variable_vec The vector to determine category count
#' @param cutoff The cutoff of the unique entries to length ratio 
#' 
#' @return Returns number of categories or 1 if continuous or gaussian
#'
get_num_categories <- function(variable_vec, cutoff=0.05){
  ifelse((length(unique(variable_vec))/length(variable_vec)) <= 0.05, length(unique(variable_vec)), 1)
  
}

#' ## Linear Model Module ####
#' 
#' ## Degrees of Freedom ##
#' #' get_dof
#' #' 
#' #' Calculates the degrees of freedom based on the phenotype vector
#' #' If the phenotype is binary, degrees of freedom is the number of samples in the
#' #' smaller group divided by 10 and the family is binomial
#' #' 
#' #' If the phenotype is continuous, degrees of freedom is n/10 and familty is "gaussian"
#' #'
#' #' @param phenotype_vec phenotype vector
#' #' 
#' #' @return Degrees of Freedom for linear model, family for linear model
#' #' 
#' get_dof <- function(
#'   phenotype_vec
#' ){
#'   
#'   # Binary phenotype
#'   if (dim(table(phenotype_vec)) == 2){
#'     dof <- round(min(table(phenotype_vec))/10)
#'     family <- "binomial"
#'   }
#'   # Continuous phenotype
#'   else { 
#'     dof <- length(phenotype_vec)/10 
#'     family <- "gaussian"
#'   }
#'   return(list(dof, family))
#' }
#' 
#' ## Lambda Calculation ##
#' #' get_lambda
#' #' 
#' #' Calculates the best value for lambda to be used in a regularized model
#' #' to maintain a certain degrees of freedom
#' #' 
#' #'
#' #' @param data_mat Matrix of data
#' #' @param dof Degrees of freedom
#' #' 
#' #' @return lambda value to be used in regularized linear model
#' #' 
#' get_lambda <- function(
#'   data_mat, 
#'   dof
#' ){
#'   d2 = svd(data_mat)$d^2
#'   ff <- function(la)  (dof-sum(sapply(d2,function(x) x/(x+la ))))^2
#'   guess= (mean(sqrt(d2))^2)/(dof/length(d2)) - mean(sqrt(d2))^2
#'   opt = optim(par = guess, ff, method = "Brent", lower = 0, upper = 1e8)
#'   
#'   structure(opt$par, d = c(dfh = sum(sapply(d2,function(x) x/(x+opt$par))), df = dof))
#' }
#' 
#' ## Module Scoring wrapper function ##
#' #' scoring_func_wrapper
#' #' 
#' #' Wrapper to calculate raw p-value of the module formed by internal node i
#' #' 
#' #'
#' #' @param i Dendrogram node index
#' #' @param R R struct
#' #' @param phenotype_vec Vector of sample phenotypes
#' #' @param confounders Vector of confounder column names
#' #' @param score_method Scoring method (either "lm" or "pc")
#' #' @param return_BIC Boolean to return BIC instead of p-value
#' #' 
#' #' @return Unadjusted p-value of module formed at node i
#' #' 
#' scoring_func_wrapper <- function(
#'   i,
#'   R,
#'   phenotype_vec,
#'   confounders,
#'   score_method,
#'   return_BIC = FALSE
#' ){
#'   # Get degrees of freedom based on number of samples
#'   # And set the family for the glm function
#'   dof <- get_dof(phenotype_vec)[[1]]
#'   family <- get_dof(phenotype_vec)[[2]]
#'   me <- R$HCL$merge
#'   
#'   # Leaf case
#'   if (i > dim(me)[1]){
#'     members <-(i - dim(me)[1])
#'     # Just return univariate model p-value or BIC
#'     return (score_regularized(R$data[,members],
#'                               phenotype_vec,
#'                               data.frame(R$samples[,confounders]),
#'                               dof,
#'                               family,
#'                               return_BIC))
#'   } 
#'   
#'   else{
#'     
#'     members <-  R$clusts[i][[1]]
#'   
#'     if (score_method == "lm"){
#'       
#'       # Number of nodes in module cannot exceed number of samples
#'       if (length(members)>length(phenotype_vec)){
#'         return(NA)
#'       }
#'       
#'       else{
#'         return (score_regularized(R$data[,members],
#'                        phenotype_vec,
#'                        data.frame(R$samples[,confounders]),
#'                        dof,
#'                        family,
#'                        return_BIC))
#'       }
#'     }
#'   
#'     if (score_method == "pc"){
#'       
#'       centered_dat <- scale(R$data[,members], center = T, scale = T)
#'       
#'       # pca
#'       pca <- prcomp(centered_dat)
#'       
#'       # Use the all the pcs up to our desired degrees of freedom
#'       pc_data<-pca$x[,1]
#'       
#'       return (score_regularized(pc_data, 
#'                               phenotype_vec, 
#'                               data.frame(R$samples[,confounders]), 
#'                               dof,
#'                               family,
#'                               return_BIC))
#'     }
#'   }
#' }
#' 
#' 
#' 
#' ## Module Scoring Function ##
#' #' score_regularized
#' #' 
#' #' Calculate raw p-value of linear model with input data and output phenotype vector
#' #' 
#' #'
#' #' @param data Data of nodes in module being tested
#' #' @param phenotype_vec Vector of sample phenotypes
#' #' @param confounder_data Data frame of confounder data
#' #' @param dof Degrees of freedom for model
#' #' @param family Family for linear model to use
#' #' @param return_BIC Boolean to return BIC instead of p-value
#' #' 
#' #' @return Unadjusted p-value of module
#' #' 
#' score_regularized <- function(
#'   data,
#'   phenotype_vec, 
#'   confounder_data,
#'   dof,
#'   family,
#'   return_BIC
#' ){
#'   
#'   #Center data
#'   centered_dat <- scale(data, center = T, scale = F)
#'   
#'   #Combine with confounders
#'   full_data<-data.frame(centered_dat, confounder_data)
#'   full_data <-as.matrix(as.data.frame(lapply(full_data, as.numeric)))
#'   
#'   # Requirement of glm
#'   if (family == "binomial"){
#'     phenotype_vec <- phenotype_vec %>% as.factor()
#'   }
#'   
#'   # Base model with only confounders (or only intercept if no confounders)
#'   if (dim(confounder_data)[2]!=0){
#'     confounder_data<- as.matrix(as.data.frame(lapply(confounder_data, as.numeric)))
#'     gn_conf <- glm(phenotype_vec~confounder_data, family =family)
#'   }
#'   else{
#'     gn_conf<-glm(phenotype_vec~1,family=family)
#'   }
#'   
#'   # Degrees of freedom exceeds features, no regularization (all pca models go here)
#'   if (dof >= ncol(full_data)){
#'     
#'     # Set degrees of freedom to that of model (with intercept) for BIC calculation
#'     dof<-ncol(full_data)+1
#'     
#'     # Model on full data 
#'     gn <- glm(phenotype_vec~full_data, family = family)
#' 
#'     # Likelihood ratio test for p-value
#'     pval <- anova(gn,gn_conf,test="LRT")$`Pr(>Chi)`[2]
#' 
#'   }
#'   
#'   # Regularization
#'   else{
#'     L <- get_lambda(centered_dat, (dof - dim(confounder_data)[2]))
#'     gn <- glmnet(x = full_data,
#'                  y = phenotype_vec,
#'                  family = family,
#'                  intercept = T,
#'                  alpha = 0,
#'                  standardize = F,
#'                  lambda = L/nrow(data),
#'                  penalty.factor = c(rep(1, dim(data)[2]), rep(0, dim(confounder_data)[2])))
#' 
#'     pval <- pchisq(deviance(gn_conf) - deviance(gn),df = (dof - (dim(confounder_data)[2])),lower.tail=F)
#' 
#'   }
#' 
#'   # BIC calculation
#'   if(return_BIC){
#'     tLL <- -deviance(gn)
#'     k <- dof
#'     n <- nobs(gn)
#'     BIC <- log(n)*k - tLL
#'     return(BIC)
#'   }
#'   pval
#' }
#'     
#'   
#' #### P-value Adjustment Module ####
#' 
#' ## P value adjustment wrapper ##
#' #' p_adjust_wrapper
#' #' 
#' #' Adjust p-values from linear models
#' #' 
#' #'
#' #' @param allpvals Pvals of all modules
#' #' @param inds Indices of allpvals that are not NA
#' #' @param R R struct
#' #' @param phenotype_vec Vector of sample phenotypes
#' #' @param confounders Vector of confounder column names
#' #' @param score_method Scoring method (either "lm" or "pc")
#' #' @param adjust_method Either "wy" for westfall-young or other adjust method accepted by p.adjust
#' #' 
#' #' @return Adjusted p-values for all modules
#' #' 
#' p_adjust_wrapper <- function(
#'   allpvals,
#'   inds,
#'   R,
#'   phenotype_vec,
#'   confounders,
#'   score_method,
#'   adjust_method,
#'   nrand
#' ){
#'   
#'   # West-fall young adjustment
#'   if (adjust_method == "wy"){
#'     return(adjust_wy(
#'       allpvals,
#'       inds,
#'       R,
#'       phenotype_vec,
#'       confounders,
#'       score_method,
#'       nrand
#'     ))
#'   }
#'   
#'   # p.adjust 
#'   else {
#'     return(p.adjust(allpvals, method = adjust_method, n = length(inds)))
#'   }
#' }
#' 
#' ## Westfall-Young pvalue adjustment ##
#' #' adjust_wy
#' #' 
#' #' Adjust p-values from linear models using Westfall-Young Randomization
#' #' 
#' #'
#' #' @param allpvals Pvals of all modules
#' #' @param inds Indices of allpvals that are not NA
#' #' @param R R struct
#' #' @param phenotype_vec Vector of sample phenotypes
#' #' @param confounders Vector of confounder column names
#' #' @param score_method Scoring method (either "lm" or "pc")
#' #' @param nrand Number of random iterations to create null p distribution
#' #' 
#' #' @return Adjusted p-values for all modules
#' #' 
#' adjust_wy <- function(
#'   allpvals,
#'   inds,
#'   R,
#'   phenotype_vec,
#'   confounders, 
#'   score_method,
#'   nrand
#' ){
#'   #### WY p-values ----
#'   ## -> randomize outcome, run all tests, record smallest p-value of each iteration
#' 
#'   wy.null <- foreach(r = 1:nrand) %dopar% {  
#'     min(sapply(inds, 
#'                scoring_func_wrapper, 
#'                R, 
#'                sample(phenotype_vec), 
#'                confounders, 
#'                score_method)) 
#'     } %>% unlist()
#'   
#'   # NAs can occur (very rarely), cut them out
#'   cut <- which(is.na(wy.null))
#'   if (length(cut)>0) {
#'     print(sprintf("NAs found in random samples: %s", paste0(cut, collapse = ',')))
#'     wy.null <- wy.null[-cut]
#'     nrand <- length(wy.null)
#'   }
#'   sapply(allpvals, function(pval){sum(wy.null<pval)}) / nrand
#' }
#' 
#' pc_plot<-function(R,phenotype_vec,num_pcs){
#'   
#'   expvars<-c()
#'   sizes<-c()
#'   for (i in 1:length(R$clusts)){
#'     members<-R$clusts[[i]]
#'     centered_dat <- scale(R$data[,members], center = T, scale = T)
#'     
#'     # pca
#'     pca <- prcomp(centered_dat)
#'     # explained variance
#'     expvar <- (pca$sdev)^2 / sum(pca$sdev^2)
#'     #list(length(members),sum(expvar[1:min(length(expvar),num_pcs)]))
#'     sizes<-c(sizes,length(members))
#'     expvars<-c(expvars,sum(expvar[1:min(length(expvar),num_pcs)]))
#'   }
#'   
#'   #sizes<-unlist(lapply(1:length(pcs),function(i) pcs[i][[1]][[1]]))
#'   #expvars<-unlist(lapply(1:length(pcs),function(i) pcs[i][[1]][[2]]))
#'   
#'   data.frame(sizes,expvars)
#' }
#' 
#' 
#' get_sig_child_density <- function(R, phenotype, confounders){
#'    sig_kids<-lapply(1:nleaves(R$HCL), function(i) {
#'     #cor.test(R$data[,i], R$samples$DIAB)$p.value
#'      data <- data.frame(pheno=R$samples[[phenotype]], molecule = R$data[,i], R$samples[,confounders])
#'      reg <- lm(pheno ~ ., data = data)
#'      summary(reg)$coefficients["molecule", "Pr(>|t|)"]
#'      
#'    }) %>% unlist %>%  p.adjust(method="bonferroni")
#' 
#'   sig_kids<-1*(sig_kids<0.05)
#'   
#'   densities <- lapply(1:5134, function(i) {
#'     members <-  R$clusts[i][[1]]
#'     sum(sig_kids[members])/length(members)
#'   }) %>% unlist()
#'   c(densities, sig_kids)
#' }
#' 
#' #dens_bins<-bin(c(densities %>% unlist, sig_kids), nbins=20)
#' 
#' #fun_color_range <- colorRampPalette(c("#1b98e0", "white"))
#' #my_colors <- fun_color_range(20)
#' 
#' 
#' get_mgm <- function(i, R, phenotype, confounders){
#'   me <- R$HCL$merge
#'   
#'   # Leaf case
#'   if (i > dim(me)[1]){
#'     return (NA)
#'   } 
#'   else{
#'     members <-  R$clusts[i][[1]]
#'     confounder_data <- R$samples[,confounders]
#'     #Combine with confounders
#'     full_data<-data.frame(R$data[,members], confounder_data)
#'     full_data <-as.matrix(as.data.frame(lapply(full_data, as.numeric)))
#'     
#'     no_nas = apply(is.na(full_data),1,sum) == 0
#'     data_no_na <- cbind(R$samples$DIAB[no_nas], full_data[no_nas,])
#'     fit_mgm <- mgm(data_no_na, type=c("c", rep("g", length(members)), "g","c","g"), level=c(2, rep(1, length(members)),1,2,1))
#'     adj_mat <- 1* (fit_mgm$pairwise$wadj!=0)
#'     nB= 100 # number of bootstrap replicates
#'     mgm_bootstrapped <- resample(fit_mgm, data_no_na, nB)
#'     
#'     binary_counts <- lapply(mgm_bootstrapped$models, function(i){
#'       1*(i$pairwise$wadj[1,]!=0)
#'     })
#'     
#'     counts_mat<-do.call(rbind,binary_counts)
#'     mgm_boots<-apply(counts_mat, 2, sum)/100
#'     
#'   }
#'   list(adj_mat, mgm_boots)
#'   
#' }
#' 
#' plot_mgm_network <- function(R, i){
#'   members <-  R$clusts[i][[1]]
#'   names <- R$annos$name[members]
#'   
#'   no_nas = apply(is.na(R$data[,members]),1,sum) == 0
#'   data_no_na <- cbind(R$samples$DIAB[no_nas], R$data[no_nas,members])
#'   fit_mgm <- mgm(data_no_na, type=c("c", rep("g", length(members))), level=c(2, rep(1, length(members))))
#'   edge_mat <- fit_mgm$pairwise$wadj
#'   adj_mat <- 1*(edge_mat!=0)
#'   adj_mat[lower.tri(adj_mat)] <- 0
#'   from_to<-which(adj_mat==1,arr.ind=T)
#'   nodes<-data.frame(id= 1:(length(members)+1), label=c("Diabetes",names))
#'   edges<-data.frame(from=from_to[,1], to=from_to[,2])
#'   
#'   rownames(adj_mat) <- c("Diabetes", members)
#'   colnames(adj_mat)<- c("Diabetes",members)
#'   G<- graph_from_adjacency_matrix(adj_mat)
#'   vs <- V(G)
#'   es <- as.data.frame(get.edgelist(G))
#'   
#'   network<-visNetwork(nodes,edges)
#'   network
#' 
#' }
#' 
