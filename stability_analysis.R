# Stability Testing

library(mgm)

# Overall model stability (same data, many runs)

# i = 2038
# members <-  QD$clusts[i][[1]]
# names <- QD$annos$name[members]
# 
# no_nas = apply(is.na(QD$data[,members]),1,sum) == 0
# data_no_na <- cbind(QD$samples$DIAB[no_nas], QD$data[no_nas,members])
# 
# data_df <- as.data.frame(data_no_na)
# names(data_df)<- c("Diabetes", names(data_df)[2:ncol(data_df)])
# data_df$Diabetes <- as.factor(data_df$Diabetes)
# 
# X <- QD$data[no_nas,members]
# Y <- QD$samples$DIAB[no_nas]
# family <- "binomial"
# 
# #MGM
# 
# library(mgm)
# mgm_runs<-lapply(1:100, function(n){
#   fit_mgm <- mgm(data_no_na, type=c("c", rep("g", length(members))), level=c(2, rep(1, length(members))))
#   fit_mgm$pairwise$wadj[1,]
# })
# 
# binary_counts <- lapply(mgm_runs, function(i){
#   1*(i!=0)
# })
# 
# counts_mat<-do.call(rbind,binary_counts)
# mgm_probs<-apply(counts_mat, 2, sum)/100
# 
# #Elastic net
# 
# library(glmnet)
# glm_runs <- lapply(1:100, function(n){
#  # Training ELastic Net Regression model
#   
#   fit = glmnet(X, Y, alpha=1, lambda=cv.glmnet(X, Y)$lambda.1se)
#   1*(c(1:length(members))%in%coef(fit)@i)
# })
# 
# counts_mat<-do.call(rbind,glm_runs)
# glm_probs <- apply(counts_mat, 2, sum)/100
# 
# probs<-data.frame(Nodes=c("Diabetes",names), GLM_probs = c("NA", glm_probs), MGM_probs = mgm_probs)
# 
# 
# 
# # Parameter estimate stability (bootstrapped data)
# 
# nB= 100 # number of bootstrap replicates
# mgm_bootstrapped <- resample(fit_mgm, data_no_na, nB)
# 
# binary_counts <- lapply(mgm_bootstrapped$models, function(i){
#   1*(i$pairwise$wadj[1,]!=0)
# })
# 
# counts_mat<-do.call(rbind,binary_counts)
# mgm_boots<-apply(counts_mat, 2, sum)/100
# 
# glm_bootstrapped <- lapply(1:nB, function(x){
#   # Training ELastic Net Regression model
#   n <- nrow(X)
#   boot_inds = sample(1:n, size = n, replace = TRUE)
#   fit = glmnet(X[boot_inds,], Y[boot_inds], alpha=1, lambda=cv.glmnet(X, Y)$lambda.1se)
#   1*(c(1:length(members))%in%coef(fit)@i)
# })
# 
# counts_mat<-do.call(rbind,glm_bootstrapped)
# glm_boots <- apply(counts_mat, 2, sum)/100
# 
# bootstrap_probs<-data.frame(Nodes=c("Diabetes",names), GLM_boostrapped = c("NA", glm_boots), MGM_bootsrapped = mgm_boots)
# 
# get_rank_plot <- function(df, title) {
#   
#   df %>% 
#     dplyr::arrange(MGM_bootsrapped) %>% 
#     dplyr::mutate(rank = row_number()) %>% 
#     ggplot() + 
#     geom_point(aes(x = rank, y = MGM_bootsrapped), size = 1.3) +
#     ggtitle(title) +
#     theme_minimal() +
#     theme(axis.text=element_text(size = 9)) +
#     theme(
#       legend.position = c(0.9, 0.6), # c(0,0) bottom left, c(1,1) top-right.
#       legend.background = element_rect(fill = "white", colour = NA)
#     )
#   
# } 

confounders_qmdiab<-c("AGE","SEX","BMI")
get_edges_linear <- function(i, R, phenotype, confounders){
  me <- R$HCL$merge
  
  # Leaf case
  if (i > dim(me)[1]){
    return (NA)
  } 
  else{
    members <-  R$clusts[i][[1]]
    confounder_data <- R$samples[,confounders]
    #Combine with confounders
    full_data<-data.frame(R$data[,members], confounder_data)
    names <- c(phenotype, names(full_data))
    full_data <-as.matrix(as.data.frame(lapply(full_data, as.numeric)))
    
    no_nas = apply(is.na(full_data),1,sum) == 0
    data_no_na <- cbind(R$samples[[phenotype]][no_nas], full_data[no_nas,])
    fit_mgm <- mgm(data_no_na, type=c("c", rep("g", length(members)), "g","c","g"), level=c(2, rep(1, length(members)),1,2,1))
    adj_mat <- 1* (fit_mgm$pairwise$wadj!=0)
    adj_mat
  }
}

get_times <- function(R){
  test_indices <- which(R$clust_info$mgm_capable[1:5134] == T)
  sizes<-c()
  times<-c()
  for (i in test_indices){
   sizes<- c(sizes, R$clust_info$Size[i]) 
   start.time <- Sys.time()
   get_edges_linear(i, R, "DIAB", confounders_qmdiab)
   end.time <- Sys.time()
   time.taken <- end.time - start.time
   times<-c(times,time.taken)
  }
  data.frame(Times = times, Sizes = sizes)
}
