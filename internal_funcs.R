### internal functions
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(shiny))


#' initialize
#' 
#' Takes as input a data matrix, annotations on the samples and annotations on the molecules
#' and creates a list structure, R with the hierarchical clustering structure and cluster informations
#'

#' @param data.matrix Matrix with samples as rows and biomolecules as columns
#' @param sample.data Matrix of annotation data of samples
#' @param mol.data Matrix of annotation data of biomolecules
#' 
#' @return R struct with the data matrix, correlation matrix, 
#' sample and molecular annotations, hierarchical structure, cluster membership
#' node order within the hierarchical structure and a color vector for the platforms

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

#' get_node_color
#' 
#' Gets the color of the nodes within the hierarchical tree based on 
#' BIC and p-value
#'
#' @param R R struct
#' @param i index of module we are coloring
#' @param signif Is the module i significant?
#' 
#' @return color of node i 
#
get_node_color <- function(
  R,
  i,
  signif,
  node_color,
  node_color_light
){
  hc <- R$HCL
  internal_nodes <- dim(hc$merge)[1]
  if (i %in% signif[!is.na(signif)]){
    if (i > (internal_nodes)) return(node_color)
    else {
      
      children <- hc$merge[i,]
      child1 <- if (children[1]<0) (internal_nodes + abs(children[1])) else children[1]
      child2 <- if (children[2]<0) (internal_nodes + abs(children[2])) else children[2]
      
      ## Children aren't significant, combination is
      if (R$pvals[child1] > 0.05 & R$pvals[child2] > 0.05) return(node_color)
      
      ## Children are significant
      if(R$pvals[child1] < 0.05 & R$pvals[child2] < 0.05){
        
        if (R$BIC[[i]] < R$BIC[[child1]] & R$BIC[[i]]<R$BIC[[child2]]) return(node_color) else return(node_color_light)
      }
      
      ## One child is significant
      else {
        sig_child <- if (R$pvals[child1] < R$pvals[child2]) child1 else child2
        
        if(R$BIC[[i]] < R$BIC[[sig_child]]) return(node_color) else return(node_color_light)
      }
    }
  }
  else{
    return("black")
  }
}

#' get_node_label
#' 
#' Gets the label of a node
#' the label is the molecule name if the node is a leaf
#' the label is the BIC and p-value if the node is internal
#'
#' @param R R struct
#' @param i index of module we are labeling
#' 
#' @return label of node i 
#' 
get_node_label <- function(
  R, 
  i
){
  internal_nodes <- dim(R$HCL$merge)[1]
  if (i > (internal_nodes)) return(R$HCL$labels[(i - internal_nodes)])
  else return(paste("BIC: ", round( R$BIC[[i]], digits = 3), ", p-value: ", R$pvals[i]))
}

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

#' get_anno_data
#' 
#' Gets the annotation data of nodes in a module
#' and structures it for a sunburst plot
#'
#' @param R R struct
#' @param i index of module we are annotating
#' 
#' @return dataframe with labels, parents, and values for sunburst plot
#' 
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

#' get_plat_colors
#' 
#' Assign a unique color to each platform in our dataset
#'
#' @param R R struct
#' @param plat_list lis of platforms in our data
#' 
#' @return vector of colors corresponding to each platform
#' 
get_plat_colors <- function(
  R,
  plat_list
){
  mapply(function(x) R$col_vector[which(R$platform == x)], plat_list)
}

#' get_plat_colors
#' 
#' Convert dendrogram to adjacency matrix for visualization
#'
#' @param hclust_obj the hclust object to convert
#' 
#' @return an adjacency matrix of the dendrogram where each parent has an edge with its children
#' 

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

#' make_R
#' 
#' Wrapper function taking in preprocessed data, outputting significantly labelled dendrogram
#'
#' @param SE SummarizedExperiment object with preprocessed data
#' @param phenotype name of the phenotype we are testing for association
#' @param confounders vector of confounder names
#' @param save_file do you want to save the R struct?
#' @param filename name of file if save_file is true
#' @param cores number of cores to use to run algorithm
#' 
#' @return R struct with significance data to be passed into shiny app
#' 

make_R <- function(
  SE,
  phenotype,
  confounders,
  node_color,
  node_color_light,
  save_file = F,
  filename = "outfile.rmd",
  cores = 6
){
  data_mat <- t(assay(SE))
  sample_data <- data.frame(colData(SE))
  mol_data <- data.frame(rowData(SE))
  R <- initialize(data_mat, sample_data, mol_data) %>% 
    find_sig_clusts(phenotype, confounders, cores, node_color, node_color_light, nrand = 50)
  if(save_file){
    save(R, file = filename)
  }
  R$phenotype = phenotype
  R
}

##################################################################

#' cluster_net
#' 
#' Create network view of module using minimum-spanning tree method
#'
#' @param R R struct
#' @param i index of parent of module to make network view
#' 
#' @return network which is a plotly network of the module
#' 

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

#' plot_dend
#' 
#' Create dendrogram view of data
#'
#' @param R R struct
#' @param order_coords coordinates of points in plot
#' 
#' @return dend_network which is a plotly network view of the dendrogram
#' 


plot_dend <- function(
  R,
  order_coords
){
  
  dend_G <- graph_from_adjacency_matrix(dend_to_adj_mat(R$HCL)) 
  es <- as.data.frame(get.edgelist(dend_G))
  
  dend_network <- 
    plot_ly(x = ~order_coords$x,
            y = ~order_coords$y) %>% 
    add_segments(data = data.frame(
      x = c(order_coords$x[es$V1], order_coords$x[es$V2]),
      xend = c(order_coords$x[es$V2],order_coords$x[es$V2]), 
      y = c(order_coords$y[es$V1], order_coords$y[es$V1]), 
      yend = c(order_coords$y[es$V1], order_coords$y[es$V2])),
      x = ~x,xend = ~xend,y = ~y,yend = ~yend,
      mode='lines',color = I('black'),size = I(1), alpha = 0.5)%>%
    add_trace(type='scatter',
              mode = "markers",
              text = ~R$labels,
              marker = list(color = "black"))#list(color = R$colors))
  dend_network
}

##################################################################



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

