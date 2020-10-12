#### This is Alex and Ben :)

require(tidyverse)
library(pheatmap)
library(Hmisc)
library(RCy3)
library(igraph)
library("dplyr")
library("psych")

mt.checkout("PACKAGE_200612")
mt.quickload()
load(data.makepath("shareddata/QMDiab/qmdiab_2019_03_13.rda"))


repeats <- function(
  D,
  i
){
  to_return <- c()
  for (j in (i+1):nrow(assay(D))){
    if(all(assay(D)[i,] == assay(D)[j,])){
      to_return <- c(to_return, j)
    }
  }
  to_return
}

remove_repeats<- function(
  D,
  mc.cores
){
  reps <- unlist(pbmclapply(1:(nrow(D)-1), function(i) repeats(D, i), mc.cores = mc.cores)) 
  D <- D[-reps,]
  D
}



preprocess_steps <- function(
  D  # SummarizedExperiment
){
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D)) 
  
  # remove samples with no data
  present = apply(is.na(assay(D)),2,sum) != nrow(D)
  D <- D[,present]
  
  # change 0s to NAs
  assay(D)[assay(D) == 0] <- NA
  
  D_alone <-

    mt_pre_filtermiss(D,metMax=0.2) %>%
    
    mt_pre_filtermiss(sampleMax=0.1) %>%
    
    mt_pre_norm_quot() %>%
    
    mt_pre_trans_log() %>%
    
    mt_pre_trans_scale() %>% 
    
    mt_pre_outliercorrection() %>% 
    
    mt_pre_filtermiss(sampleMax=0.1) %>% 
    
    mt_pre_impute_knn() 
    
  # return
  D_alone
  
}

IgA_prots <- c("ENI1", "HYT1", "IIV1", "LAGC", "LAGY", "LSL", "SES", "TPL")
preprocess_IgA <- function(
  D
){
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D)) 
  
  # remove samples with no data
  present = apply(is.na(assay(D)),2,sum) != nrow(D)
  D <- D[,present]
  
  # change 0s to NAs
  assay(D)[assay(D) == 0] <- NA
  
  D_alone <-
    mt_pre_filtermiss(D, metMax=0.2) %>%
    
    mt_pre_filtermiss(sampleMax=0.1) 
  for (prot in IgA_prots){
    indices <- which(startsWith(rowData(D_alone)$name, prot))
    if (length(indices)!= 0){
      temp_D <- D_alone[indices[1]:indices[length(indices)],]
      if (exists("dfTot")){
        dfTot <-  rbind(dfTot, mt_pre_norm_quot(temp_D))
      }
      else{
        dfTot <- mt_pre_norm_quot(temp_D)
        
      }
    }
  }
  
  D_alone <-
    dfTot %>%
    
    mt_pre_trans_log() %>%
    
    mt_pre_impute_knn()
  # return
  D_alone
}



bind_SE_no_NA <- function(
  platforms  # SummarizedExperiments to combine
){
  for (platform in platforms){
    if (exists("dfTot")){
      missing1 <- setdiff(colnames(rowData(dfTot)), colnames(rowData(platform)))
      missing2 <- setdiff(colnames(rowData(platform)), colnames(rowData(dfTot)))
      for (colname in missing1){
        rowData(platform)[[colname]] <- rep(NA, nrow(platform))
      }
      for (colname in missing2){
        rowData(dfTot)[[colname]]<- rep(NA, nrow(dfTot))
      }
      
      overlap <- intersect(colnames(dfTot), colnames(platform))
      dfTot <-  rbind(dfTot[,overlap], platform[,overlap])
    }
    else{
      dfTot <-platform
    }
  }
  dfTot
}



find_missing_samples <- function(
  datasets  # List of datasets in qmdiab
){
  output <- matrix(ncol=length(datasets), nrow=2 )
  colnames(output) <- datasets
  rownames(output) <- c("Present", "Missing")
  heatmapvec <- matrix(nrow = ncol(qmdiab[[datasets[1]]]), ncol = length(datasets))
  colnames(heatmapvec)<- datasets
  
  for (i in 1:length(datasets)){
    
    D <- qmdiab[[datasets[i]]]
    present = apply(is.na(assay(D)),2,sum) != nrow(D)
    output[1,i] <- count(present)
    output[2,i] <- ncol(D) - count(present)
    heatmapvec[,i] <- as.vector(present)
  }
  
  heatmapvec <- apply(heatmapvec, c(1,2), as.numeric)
  heatmap(heatmapvec, Rowv=NA, Colv=NA, col = c("Black", "White"), xlab = "Platform", ylab = "sampleID", main = "Missing Samples" )
  
}


find_overlap <- function(
  datasets     # List of datasets in qmdiab
){
  samples_list <- c()
  metabolites_list <- c()
  datasets_list <- list()
  #for (i in 1:length(datasets)){
  combos <- combn(datasets, 5)
  print(combos)
  for (i in 1:ncol(combos)){
    
    bound <- bind_platforms(as.vector(combos[,i]))
    bound_no_NA <- bound[ , colSums(is.na(bound)) == 0]
    
    samples_list <- c(samples_list, ncol(bound_no_NA))
    metabolites_list <- c(metabolites_list, nrow(bound_no_NA))
    datasets_list <- append(datasets_list, list(combos[,i]))
    print(datasets_list)
    print (paste(i, " out of ", ncol(combos), " have run"))
    #}
  }
  final <- list(samples_list, metabolites_list, datasets_list)
  final
}


plotToCytoscape <- function(
  D,
  matrix, # Correlation, p-value, or adjacency matrix
  cutoff = NA, # Cutoff for the values, NA if the matrix is an adjacency matrix
  less_than = F # if the values of the matrix have to be less than the cutoff
){
  significant <- matrix
  if (!is.na(cutoff)){
    if (less_than){
      significant <-matrix <= cutoff
    }
    else{
      significant <-matrix >= cutoff
    }
  }
  
  significant_numeric <- 1*significant
  significant_numeric[lower.tri(significant_numeric)] <- 0
  net <- igraph::graph_from_adjacency_matrix(significant_numeric)
  
  V(net)$name <- rowData(D)$name
  V(net)$platform <- rowData(D)$platform
  net
  #createNetworkFromIgraph(net)
  #significant_numeric
}

process_PM = preprocess_steps(qmdiab$PM)
process_UM = preprocess_steps(qmdiab$UM)
process_SM = preprocess_steps(qmdiab$SM)
process_HDF = preprocess_steps(qmdiab$HDF)
process_CM = preprocess_steps(qmdiab$CM)
process_BM = preprocess_steps(qmdiab$BM)
process_BRAIN = preprocess_steps(qmdiab$BRAIN)
process_SOMA = preprocess_steps(qmdiab$SOMA)
process_LD = preprocess_steps(qmdiab$LD)
process_IgA = preprocess_IgA(qmdiab$IgA)
process_GP = preprocess_steps(qmdiab$rawGP)
process_IgG = preprocess_steps(qmdiab$rawIgG)


platform_list <- list(process_PM, process_UM, process_SM, process_HDF, process_CM, process_BM, process_BRAIN, process_SOMA, process_LD, process_IgA, process_GP, process_IgG)
all_platforms <- bind_SE_no_NA(platform_list)
#ks <- lapply(1:nrow(all_platforms), function(x) ks.test(scale(assay(all_platforms)[x,]), pnorm)$p.value)
#names <- c("PM", "UM", "SM", "HDF", "CM","BM","BRAIN","SOMA","LD", "IgA","GP","IgG")
# samples available in all omics 
#sample_ids = lapply(all_platforms, colnames) %>% unlist %>% table %>% `==`(length(all_platforms)) %>% which %>% names
# keep those samples available in all omics
# nods = lapply(all_platforms, function(x) x[,sample_ids])
# par(mfrow = c(4,6), mar=c(1,1,1,1))
# for (i in 1:length(nods)){
#   print (i)
#   pc1 <- prcomp(t(assay(nods[[i]])) %>% apply(2, function(a){a[is.na(a)]=min(a,na.rm = T);a}))$x[,"PC1"] 
#   pc.ks <- ks.test(scale(pc1), pnorm)$p.value
#   hist(pc1, 20,  main = paste(names[i], "pval: ", round(pc.ks, digits = 2)),freq = F, col = "gray")
#   lines(density(pc1), col  = "red")
#   first.ks <- ks.test(scale(as.vector(assay(nods[[i]][prcomp( t(assay(nods[[i]])) )$rotation[,"PC1"] %>% abs %>% which.max,]))), pnorm)$p.value
#   assay(nods[[i]][prcomp( t(assay(nods[[i]])) )$rotation[,"PC1"] %>% abs %>% which.max,]) %>% {
#     hist(.,20,freq = F, main = paste("Main contributor to", names[i], "pval: ", round(first.ks, digits = 2)), col="gray");
#     density(.) %>% lines(col = "red")
#   }
# }

#metab <- bind_SE_no_NA(list(process_HDF, process_UM, process_SM))

### ADNI DATA LOAD ###
# lipids <-read.csv(data.makepath("ADNI/adni1.lipids.medadjusted.csv"))
# meikle <-read.csv(data.makepath("ADNI/adni1.meikle.medadjusted.csv"))
# p180 <-read.csv(data.makepath("ADNI/adni1go2.p180.medadjusted.csv")) %>% dplyr::select(-Cohort) 
# ba <-read.csv(data.makepath("ADNI/adni1go2.ba.medadjusted.csv")) %>% dplyr::select(-Cohort) 
# covars <- read.csv(data.makepath('ADNI/adni1go2.phenotypes.covariates.csv'))
# 
# 
# adni_join <- inner_join(lipids, meikle, p180, ba, by = "RID")
# platform <- c(rep("lipids", ncol(lipids)-1), rep("meikle", ncol(meikle)-1), rep("p180", ncol(p180)-1), rep("ba", ncol(ba)-1))
# 
# multi_inner <- Reduce(
#   function(x, y, ...) merge(x, y, ...), 
#   list(lipids, meikle, p180, ba)
# )
# 
# ADNI_SE <- SummarizedExperiment(assays = t(multi_inner %>% dplyr::select(-RID)), 
#                                 colData = covars %>% filter(RID %in% multi_inner$RID), 
#                                 rowData = data.frame(
#                                   name = colnames(multi_inner)[2:length(colnames(multi_inner))], platform))
# 

# locs <- unlist(strsplit(rowData(all_platforms)$HMDB, ","))
# locs <- unique(unlist(strsplit(locs, "\\|")))
# locs <- locs[!is.na(locs)]
# locs <- locs[locs!= 'NA']
# 
# loc_freq <- c()
# for (loc in locs){
#   inds <- grep(loc, rowData(all_platforms)$HMDB)
#   loc_freq<- c(loc_freq, length(unique(rowData(all_platforms)$platform[inds])))
# }
# geneNet <- geneNetMatrix(t(assay(all_platforms)))
# pcors <- geneNet$pcor
# pvals <- geneNet$pval
# diag(pvals)<-1
# within <- c()
# platform_names <- unique(rowData(all_platforms)$platform)
# adj_p <- p.adjust(as.vector(pvals[upper.tri(pvals)]), method = "bonferroni", n = length(pvals[upper.tri(pvals)]))
# pvals_mat <- matrix(1, nrow = nrow(pvals), ncol = ncol(pvals))
# pvals_mat[upper.tri(pvals_mat)]<-adj_p
# pvals_mat[lower.tri(pvals_mat)] <- adj_p
# cutoff <- -log10(pvals[which(pvals_mat == max(pvals_mat[pvals_mat < 0.05]), arr.ind = T)])
# mycor <- rcorr(t(assay(all_platforms)), type="pearson")
# diag(mycor$P)<-1
# averages <- c()
# for (plat in platform_names){
#   for (plat2 in platform_names){
#     inds1 <- which(rowData(all_platforms)$platform == plat)
#     inds2 <- which(rowData(all_platforms)$platform == plat2)
#     add <- as.vector(mycor$P[inds1, inds2])
#     averages <- c(averages, mean(-log10(add)))
#   }
# }
# plat_hm <- matrix(averages, nrow = length(platform_names), ncol = length(platform_names))
# rownames(plat_hm) <- platform_names
# colnames(plat_hm) <- platform_names
# melted <- melt(plat_hm)
# p <- ggplot(melted, aes(Var1, Var2)) + geom_tile(aes(fill = value),
#          colour = "white") + scale_fill_gradient(low = "white",high = "steelblue")
# #between <- as.vector(pvals[pvals!=0])
# 
# #within <- data.frame(edges = -log10(within))
# #between <- data.frame(edges = -log10(between))
# 
# #within$label <- "within"
# #between$label <- "between"
# 
# 
# 
# 
# 
# edge_df <- rbind(within, between)
# ggplot(edge_df[edge_df$edges > 8,], aes(edges, fill = label))+ geom_density(alpha = 0.2) + geom_vline(xintercept = cutoff)
# 
# 
# 
# 
