### Frontend Modules ###


suppressPackageStartupMessages(library(RColorBrewer))


#### Make Dendrogram ####

#### Label nodes ####
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


#' plot_dend
#' 
#' Create dendrogram view of data
#'
#' @param R R struct
#' @param order_coords coordinates of points in plot
#' @param i index of selected node, default NA
#' 
#' @return dend_network which is a plotly network view of the dendrogram
#' 


plot_dend <- function(
  R,
  order_coords,
  i
){
  
  dend_G <- graph_from_adjacency_matrix(dend_to_adj_mat(R$HCL)) 
  es <- as.data.frame(get.edgelist(dend_G))
  colors <- R$colors
  colors[i] <- "yellow"
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
              text = mapply(function(i) get_node_label(R, i), 1:nnodes(R$HCL)),
              marker = list(color = colors, size = 9))
  dend_network
}


plot_dend_fast<-function(
  R,
  order_coords,
  dend_data
){
  colors <- R$colors
    plot_ly(x = ~order_coords$x,
            y = ~order_coords$y) %>% 
    add_segments(
      x=dend_data$segments$x, 
      xend=dend_data$segments$xend, 
      y=dend_data$segments$y, 
      yend=dend_data$segments$yend)%>%
    add_trace(type='scatter',
              mode = "markers",
              text = mapply(function(i) get_node_label(R, i), 1:nnodes(R$HCL)),
              marker = list(color = colors, size = 9))
}
#### Make network view ####

#### Get platform colors ####
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
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, 
                               qual_col_pals$maxcolors, 
                               rownames(qual_col_pals)))[1:length(R$platforms)]
  mapply(function(x) col_vector[which(R$platform == x)], plat_list)
}

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
  L <- layout_with_gem(G)
  
  network <- plot_ly(x = ~L[,1], y = ~L[,2]) %>% 
    
    add_segments(data = data.frame(x = L[,1][es$V1], 
                                   xend = L[,1][es$V2], 
                                   y = L[,2][es$V1], 
                                   yend = L[,2][es$V2]),
                 x = ~x, xend = ~xend, y = ~y, yend = ~yend,
                 mode='lines',
                 color = I('black'),
                 size = I(1), 
                 alpha = 0.5) %>%
    
    add_trace(type = "scatter",
              mode = "markers", 
              text = names, 
              hoverinfo = "text",
              marker = list(color = get_plat_colors(R, platforms), 
                            size = 40))
  network
}

#### Convert dendrogram to adjacency matrix ####
#' dend_to_adj_mat
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


#### Make list view ####

### Make annotation view ###

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
  mat <- subset( R$annos[R$clusts[[i]],], select = c(platform, SUPER_PATHWAY, SUB_PATHWAY))
  labels <- c()
  parents <- c()
  values <- c()
  for (i in 1:length(colnames(mat))){
    categories <- unique(mat[,i][!is.na(mat[,i])])
    if (length(categories) != 0){
      labels <- c(labels, colnames(mat)[i], categories)
      parents <- c(parents, "", rep(colnames(mat)[i], times = length(categories)))
      values <- c(values, 0, table(mat[,i])[categories])
    }
  }
  sun_df <- data.frame(labels, parents, values)
  sun_df
}
