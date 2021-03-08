## Backend Modules ###
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(glmnet))
library(RhpcBLASctl)
blas_set_num_threads(1)
## Hierarchical Clustering Module ####

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

## Get dendrogram from best clustering method ##
#' get_cluster_method
#' 
#' Searches through all hierarchical clustering methods for one that
#' results in a dendrogram with highest cophenetic correlation to original data mat
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

## Linear Model Module ####

## Degrees of Freedom ##
#' get_dof
#' 
#' Calculates the degrees of freedom based on the phenotype vector
#' If the phenotype is binary, degrees of freedom is the number of samples in the
#' smaller group divided by 10 and the family is binomial
#' 
#' If the phenotype is continuous, degrees of freedom is n/10 and familty is "gaussian"
#'
#' @param phenotype_vec phenotype vector
#' 
#' @return Degrees of Freedom for linear model, family for linear model
#' 
get_dof <- function(
  phenotype_vec
){
  
  # Binary phenotype
  if (dim(table(phenotype_vec)) == 2){
    dof <- round(min(table(phenotype_vec))/10)
    family <- "binomial"
  }
  # Continuous phenotype
  else { 
    dof <- length(phenotype_vec)/10 
    family <- "gaussian"
  }
  return(list(dof, family))
}

## Lambda Calculation ##
#' get_lambda
#' 
#' Calculates the best value for lambda to be used in a regularized model
#' to maintain a certain degrees of freedom
#' 
#'
#' @param data_mat Matrix of data
#' @param dof Degrees of freedom
#' 
#' @return lambda value to be used in regularized linear model
#' 
get_lambda <- function(
  data_mat, 
  dof
){
  d2 = svd(data_mat)$d^2
  ff <- function(la)  (dof-sum(sapply(d2,function(x) x/(x+la ))))^2
  guess= (mean(sqrt(d2))^2)/(dof/length(d2)) - mean(sqrt(d2))^2
  opt = optim(par = guess, ff, method = "Brent", lower = 0, upper = 1e8)
  
  structure(opt$par, d = c(dfh = sum(sapply(d2,function(x) x/(x+opt$par))), df = dof))
}

## Module Scoring wrapper function ##
#' scoring_func_wrapper
#' 
#' Wrapper to calculate raw p-value of the module formed by internal node i
#' 
#'
#' @param i Dendrogram node index
#' @param R R struct
#' @param phenotype_vec Vector of sample phenotypes
#' @param confounders Vector of confounder column names
#' @param score_method Scoring method (either "lm" or "pc")
#' @param return_BIC Boolean to return BIC instead of p-value
#' 
#' @return Unadjusted p-value of module formed at node i
#' 
scoring_func_wrapper <- function(
  i,
  R,
  phenotype_vec,
  confounders,
  score_method,
  return_BIC = FALSE
){
  # Get degrees of freedom based on number of samples
  # And set the family for the glm function
  dof <- get_dof(phenotype_vec)[[1]]
  family <- get_dof(phenotype_vec)[[2]]
  me <- R$HCL$merge
  
  # Leaf case
  if (i > dim(me)[1]){
    members <-(i - dim(me)[1])
  } 
  
  else{
    
    members <-  R$clusts[i][[1]]
  }
  
  if (score_method == "lm"){
    
    # Number of nodes in module cannot exceed number of samples
    if (length(members)>length(phenotype_vec)){
      return(NA)
    }
    
    else{
      return (score_regularized(R$data[,members],
                     phenotype_vec,
                     data.frame(R$samples[,confounders]),
                     dof,
                     family,
                     return_BIC))
    }
  }

  if (score_method == "pc"){
    
    # desired explained variance
    dexpvar <- 0.95
    
    # pca
    pca <- prcomp(R$data[,members])
    # explained variance
    expvar <- (pca$sdev)^2 / sum(pca$sdev^2)
    
    # find number of components when desired explained variance is reached
    reachp <- which(cumsum(expvar)>=dexpvar)[1]
    
    pc_data <- pca$x[,reachp]
    
    # Size of module less than degrees of freedom, no regularization
   return (score_regularized(pc_data, 
                            phenotype_vec, 
                            as.matrix(R$samples[,confounders]), 
                            dof,
                            family,
                            return_BIC))
  }
}
## Module Scoring Function ##
#' score_regularized
#' 
#' Calculate raw p-value of linear model with input data and output phenotype vector
#' 
#'
#' @param data Data of nodes in module being tested
#' @param phenotype_vec Vector of sample phenotypes
#' @param confounders Vector of confounder column names
#' @param dof Degrees of freedom for model
#' @param family Family for linear model to use
#' @param return_BIC Boolean to return BIC instead of p-value
#' 
#' @return Unadjusted p-value of module
#' 
score_regularized <- function(
  data,
  phenotype_vec, 
  confounders,
  dof,
  family,
  return_BIC
){
  
  #Center data
  centered_dat <- scale(data, center = T, scale = F)
  full_data<-data.frame(centered_dat, confounders)
  full_data <-as.matrix(as.data.frame(lapply(full_data, as.numeric)))
  if (family == "binomial"){
    phenotype_vec <- phenotype_vec %>% as.factor()
  }
  
  no_na<-complete.cases(full_data,phenotype_vec)
  full_data<-full_data[no_na,] %>% as.matrix()
  confounders<-confounders[no_na,]
  phenotype_vec<-phenotype_vec[no_na]
  
  # Base model with only confounders
  if (dim(confounders)[2]!=0){
    confounders<- as.matrix(as.data.frame(lapply(confounders, as.numeric)))
    gn_conf <- glm(phenotype_vec~confounders-1, family =family)
  }
  else{
    gn_conf<-glm(phenotype_vec~-1,family=family)
  }
  
  # Degrees of freedom exceeds features, no regularization
  if (dof >= ncol(full_data)){
    gn <- glm(phenotype_vec~full_data-1, family = family)

    pval <- anova(gn,gn_conf,test="LRT")$`Pr(>Chi)`[2]
    #pval <- -expm1(pchisq(deviance(gn) - deviance(gn_conf),df = (dof - dim(confounders)[2]),log.p = T))
    
  }
  
  # Regularization
  else{
    
    L <- get_lambda(centered_dat, (dof - dim(confounders)[2]))
    gn <- glmnet(x = full_data, 
                 y = phenotype_vec, 
                 family = family, 
                 intercept = F,
                 alpha = 0, 
                 standardize = F,
                 lambda = L/nrow(data),
                 penalty.factor = c(rep(1, dim(data)[2]), rep(0, dim(confounders)[2])))

    pval <- -expm1(pchisq(deviance(gn) - deviance(gn_conf),df = (dof - dim(confounders)[2]),log.p = T))

  }
  
  # BIC calculation
  if(return_BIC){
    tLL <- -deviance(gn)
    k <- dof
    n <- nobs(gn)
    BIC<-log(n)*k - tLL
    return(BIC)
  }
  pval
}
    
  
#### P-value Adjustment Module ####

## P value adjustment wrapper ##
#' p_adjust_wrapper
#' 
#' Adjust p-values from linear models
#' 
#'
#' @param allpvals Pvals of all modules
#' @param inds Indices of allpvals that are not NA
#' @param R R struct
#' @param phenotype_vec Vector of sample phenotypes
#' @param confounders Vector of confounder column names
#' @param score_method Scoring method (either "lm" or "pc")
#' @param adjust_method Either "wy" for westfall-young or other adjust method accepted by p.adjust
#' 
#' @return Adjusted p-values for all modules
#' 
p_adjust_wrapper <- function(
  allpvals,
  inds,
  R,
  phenotype_vec,
  confounders,
  score_method,
  adjust_method
){
  
  # West-fall young adjustment
  if (adjust_method == "wy"){
    return(adjust_wy(
      allpvals,
      inds,
      R,
      phenotype_vec,
      confounders,
      score_method
    ))
  }
  
  # p.adjust 
  else {
    return(p.adjust(allpvals, method = adjust_method, n = length(inds)))
  }
}

## Westfall-Young pvalue adjustment ##
#' adjust_wy
#' 
#' Adjust p-values from linear models using Westfall-Young Randomization
#' 
#'
#' @param allpvals Pvals of all modules
#' @param inds Indices of allpvals that are not NA
#' @param R R struct
#' @param phenotype_vec Vector of sample phenotypes
#' @param confounders Vector of confounder column names
#' @param score_method Scoring method (either "lm" or "pc")
#' @param nrand Number of random iterations to create null p distribution
#' 
#' @return Adjusted p-values for all modules
#' 
adjust_wy <- function(
  allpvals,
  inds,
  R,
  phenotype_vec,
  confounders, 
  score_method,
  nrand = 1000
){
  #### WY p-values ----
  ## -> randomize outcome, run all tests, record smallest p-value of each iteration

  wy.null <- foreach(r = 1:nrand) %dopar% {  
    min(sapply(inds, 
               scoring_func_wrapper, 
               R, 
               sample(phenotype_vec), 
               confounders, 
               score_method)) 
    } %>% unlist()
  
  # NAs can occur (very rarely), cut them out
  cut <- which(is.na(wy.null))
  if (length(cut)>0) {
    print(sprintf("NAs found in random samples: %s", paste0(cut, collapse = ',')))
    wy.null <- wy.null[-cut]
    nrand <- length(wy.null)
  }
  sapply(allpvals, function(pval){sum(wy.null<pval)}) / nrand
}
