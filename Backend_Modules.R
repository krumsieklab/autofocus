## Backend Modules ###
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


get_dendro <- function(
  R,
  method = NA,
  cores = cores
){
  dist_mat <- abs_cor_dist(R)
  if (!(is.na(method))){
    return(hclust(dist_mat, method = method))
  }
  else {
    return(get_cluster_method(dist_mat, cores = cores))
  }
}


grid_search_clust <- function(
  dist_mat, 
  cores = 1
){
  
  test_vals <- seq(0,1,by=0.05)
  flex_known_cors <- unlist(mclapply(test_vals, function(i){
    hc = as.hclust(cluster::agnes(dist_mat, diss = T,method = "flexible", par.method = i))
    cor(as.vector(as.dist(cophenetic(hc))), as.vector(dist_mat), method = "spearman")
  }, mc.cores = cores))
  flex_known_data <- data.frame(flex_known_cors, test_vals,clustering_type = rep("flexible", length(flex_known_cors)))
  
  
  
  gav_known_cors <- unlist(mclapply(test_vals, function(i){
    hc = as.hclust(cluster::agnes(dist_mat, diss = T,method = "gaverage", par.method = i))
    cor(as.vector(as.dist(cophenetic(hc))), as.vector(dist_mat), method = "spearman")
  }, mc.cores = cores))
  gav_known_data <- data.frame(gav_known_cors, test_vals, clustering_type = rep("gaverage", length(gav_known_cors)))
  
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
                     as.matrix(R$samples[,confounders]),
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
  
score_regularized <- function(
  data,
  phenotype_vec, 
  confounders,
  dof,
  family,
  return_BIC
){
  
  centered_dat <- scale(data, center = T, scale = F)
  
  full_data <- cbind(centered_dat, confounders) %>% as.matrix()
  
  if (family == "binomial"){
    phenotype_vec <- phenotype_vec %>% as.factor()
  }
  
  if (dof >= ncol(full_data)){
    gn <- glmnet(x = full_data, 
               y = phenotype_vec,
               family = family, 
               intercept = F,
               alpha = 0,
               standardize = F,
               lambda = 0)
  }
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
  }
  gn_conf <- glmnet(x = confounders,
                    y = phenotype_vec,
                    family =family,
                    intercept = F,
                    lambda = 0,
                    standardize = F)
  
  if(return_BIC){
    tLL <- gn$nulldev - deviance(gn)
    k <- (dof - dim(confounders)[2])
    n <- gn$nobs
    BIC<-log(n)*k - tLL
    return(BIC)
  }
  
  # p-val of model differences
  d_dev = deviance(gn_conf) - deviance(gn)
  print(-expm1( pchisq( d_dev,df = (dof - dim(confounders)[2]), log.p = T)))
  -expm1( pchisq( d_dev,df = (dof - dim(confounders)[2]), log.p = T))
}
    
  
#### P-value Adjustment Module ####

p_adjust_wrapper <- function(
  allpvals,
  inds,
  R,
  phenotype_vec,
  confounders,
  score_method,
  adjust_method
){
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
  else {
    return(p.adjust(allpvals, method = adjust_method, n = length(inds)))
  }
}

adjust_wy <- function(
  allpvals,
  inds,
  R,
  phenotype_vec,
  confounders, 
  method,
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
               method)) 
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
