### internal functions
library(dendextend)
library(magrittr)
library(parallel)
library(foreach)
library(RColorBrewer)
library(plotly)


initialize <- function(
  data.matrix, # Matrix of raw data (samples * nodes)
  sample.data, # Matrix of sample annotations
  mol.data # Matrix of node annotations
){
  
  R <- list(data = data.matrix, samples = sample.data, annos = mol.data)
  R$C <- cor(data.matrix)
  
  R$HCL <- R %>% abs_cor_dist() %>% hclust(method = "average")
  R$order <- get_dend_indices(R, dim(R$HCL$merge)[1], c())
  
  
  R$clusts <- lapply(1:dim(R$HCL$merge)[1], function(x) get_members(R, x))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  R$col_vector = unlist(mapply(brewer.pal, 
                               qual_col_pals$maxcolors, 
                               rownames(qual_col_pals)))[1:length(R$platforms)]
  R
}

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

abs_cor_dist <- function(
  R
){
  as.dist((1-abs(R$C))) # abs correlation
}


get_anno_data <- function(
  R,
  i
){
  mat <- subset( R$annos, select = -name )
  labels <- c()
  parents <- c()
  values <- c()
  for (i in 1:length(colnames(mat))){
    categories <- unique(mat[,i][!is.na(mat[,i])])
    labels <- c(labels, colnames(mat)[i], categories)
    parents <- c(parents, "", rep(colnames(mat)[i], times = length(categories)))
    values <- c(values, 0, table(mat[,i])[categories])
  }
  sun_df <- data.frame(labels, parents, values)
  sun_df
}


get_plat_colors <- function(
  R,
  plat_list
){
  mapply(function(x) R$col_vector[which(R$platform == x)], plat_list)
}


dend_to_adj_mat <- function(
  hclust_obj
){
  
  adj_mat <- matrix(0L, nrow = (2*dim(hclust_obj$merge)[1])+1, ncol = (2*dim(hclust_obj$merge)[1])+1)
  for(i in rev(1:dim(hclust_obj$merge)[1])){
    parent <- i
    right_index <- hclust_obj$merge[i,1]
    if (right_index < 0){
      right_index <- abs(right_index) + dim(hclust_obj$merge)[1]
    }
    left_index <- hclust_obj$merge[i,2]
    if (left_index < 0){
      left_index <- abs(left_index) + dim(hclust_obj$merge)[1]
    }
    adj_mat[i, right_index] <- 1
    adj_mat[i, left_index] <- 1
  }
  adj_mat
}


get_dend_indices <- function(R,number, indices){
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


cluster_net <- function(
  R,
  i
) {
  
  members <- R$clusts[[i]]
  names <- R$annos$name[members]
  platforms <- R$annos$platform[members]
  cor_vals <- R$C[members,members]
  full_conn <- graph_from_adjacency_matrix((-1*abs(cor_vals)), 
                                           mode='undirected', 
                                           weighted = T)
  min_span <- mst(full_conn, weights = E(full_conn)$weight)
  cutoff <- abs(max(E(min_span)$weight))
  adj_mat <- 1*(abs(cor_vals) >= cutoff)
  adj_mat[lower.tri(adj_mat)] <- 0
  rownames(adj_mat) <- members
  colnames(adj_mat)<- members
  G<- graph_from_adjacency_matrix(adj_mat)
  vs <- V(G)
  vs$platform <- platforms
  es <- as.data.frame(get.edgelist(G))
  L <- layout_nicely(G)

  network <- plot_ly(x = ~L[,1], y = ~L[,2]) %>% 
  
              add_segments(data = data.frame(x = L[,1][es$V1], 
                                             xend = L[,1][es$V2], 
                                             y = L[,2][es$V1], 
                                             yend = L[,2][es$V2]),
                            x = ~x, xend = ~xend, y = ~y, yend = ~yend,
                           mode='lines',color = I('black'),size = I(1), alpha = 0.5) %>% 
              add_trace(type = "scatter",
                        mode = "markers", text = names, hoverinfo = "text",
                        marker = list(get_plat_colors(R, platforms)))
  network
}






get_cluster_method <- function(
  R
){
  known_cors <- c()
  methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median",'centroid')
  for (method in methods){
    clust <- hclust(mydist(R), method = method)
    known_cors <- c(known_cors, cor(cophenetic(clust), mydist(R), method = "spearman"))
  }
  
  known_data <- data.frame(known_cors, methods)
  flex_weighted <-c()
  flex_unweighted <-c()
  
  x <- seq(0,1,by=0.1)
  for (i in x){
    clus <- cluster::agnes(mydist(R), method = 'flexible', par.method = i)
    flex_weighted <- c(flex_weighted, cor(cophenetic(clus), mydist(X), method = "spearman"))
    
    clus <- cluster::agnes(mydist(R), method = 'gaverage', par.method = i)
    flex_unweighted <- c(flex_unweighted, cor(cophenetic(clus), mydist(R), method = "spearman"))
  }
  flex_data <- data.frame(x, flex_weighted, flex_unweighted)
  print(ggplot(data=flex_data) + 
          geom_point(aes(x=x, y=flex_weighted,color = "Weighted")) +
          geom_point(aes(x=x, y=flex_unweighted, color = "Unweighted")) +
          geom_hline(data = known_data, aes(yintercept=known_cors, color = methods))+
          ggtitle("Correlation between linkage method and distance matrix") +
          xlab("x (flexible parameter)") + ylab("Spearman correlation"))
  
}


dend_to_text_wrapper <- function(
  R,
  i,
  filename,
  list_children = F
){
  fileConn<-file(filename)
  writeLines(paste("<root>",dend_to_text(R, i, list_children),"</root>"), fileConn)
  close(fileConn)
}

dend_to_text <- function(
  R,
  i,
  list_children
){
  me <- R$HCL$merge
  #Leaf case
  if(i <0){
    index <- abs(i) + dim(me)[1]
    lines <- paste("<Name>", R$annos$name[abs(i)], "</Name>", 
               "<Type> Leaf </Type>", 
               "<pval>", R$pvals[index], "</pval>",
               "<BIC>", round(R$BIC[index], digits = 6), "</BIC>")
  }
  
  else{
    if (list_children){
      child_add <- paste("<All_children>",paste(R$annos$name[abs(get_members(R, i))], collapse= ", "),'</All_children>')
    }
    else{ child_add = ""}
    left_text <- dend_to_text(R, me[i,][1])
    right_text <-dend_to_text(R, me[i,][2])
    lines <- paste("<Name>", paste(get_anno_name(R,i),  collapse =","), "</Name>", 
               "<Type> Internal Node </Type>", 
               "<pval>", R$pvals[i], "</pval>",
               "<BIC>", round(R$BIC[i], digits = 6), "</BIC>",
               child_add,
               "<Child>", left_text, "</Child>",
               "<Child>", right_text, "</Child>")
  }
  lines
}


phenotype_correlations <- function(
  R,
  phen_start,
  phen_end
){
  phens_df <-  R$samples[, phen_start:phen_end]
  cor_mat <- matrix(0L, nrow = ncol(phens_df), ncol = ncol(phens_df))
  for (i in 1:ncol(phens_df)){
    for (j in 1:ncol(phens_df)){
      overlap <- intersect(which(!is.na(phens_df[,i])), which(!is.na(phens_df[,j])))
      cor_mat[(i-5), (j-5)]<- cor(phens_df[overlap,i],phens_df[overlap,j])
    }
  }
  rownames(cor_mat) <- colnames(R$samples)[phen_start:phen_end]
  colnames(cor_mat)<-colnames(R$samples)[phen_start:phen_end]
  cor_mat
}

