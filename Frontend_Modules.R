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
  else return(paste("Cluster ID:",R$clust_info$ClusterID[[i]] ,"BIC: ", round( R$clust_info$BIC[[i]], digits = 3), ", p-value: ", R$clust_info$pvals[i]))
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


plot_dend<-function(
  R
){
  R$clust_info %>% plot_ly(x = ~Coord_X,
          y = ~Coord_Y) %>% 
  add_segments(
    x=R$dend_data$segments$x, 
    xend=R$dend_data$segments$xend, 
    y=R$dend_data$segments$y, 
    yend=R$dend_data$segments$yend)%>%
  add_trace(type='scatter',
            mode = "markers",
            text = mapply(function(i) get_node_label(R, i), 1:nnodes(R$HCL)),
            marker = list(color = ~colors, size = 9))
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
  tic()
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
  diag(adj_mat)<-0
  from_to<-which(adj_mat==1,arr.ind=T)
  nodes<-data.frame(id=1:length(members), label=names)
  edges<-data.frame(from=from_to[,1], to=from_to[,2])
  toc()
  # tic()
  # rownames(adj_mat) <- members
  # colnames(adj_mat)<- members
  # G<- graph_from_adjacency_matrix(adj_mat)
  # vs <- V(G)
  # vs$platform <- platforms
  # es <- as.data.frame(get.edgelist(G))
  # toc()
  tic()
  network<-visNetwork(nodes,edges)
  # L <- layout_with_gem(G)
  # 
  # network <- plot_ly(x = ~L[,1], y = ~L[,2]) %>% 
  #   
  #   add_segments(data = data.frame(x = L[,1][es$V1], 
  #                                  xend = L[,1][es$V2], 
  #                                  y = L[,2][es$V1], 
  #                                  yend = L[,2][es$V2]),
  #                x = ~x, xend = ~xend, y = ~y, yend = ~yend,
  #                mode='lines',
  #                color = I('black'),
  #                size = I(1), 
  #                alpha = 0.5) %>%
  #   
  #   add_trace(type = "scatter",
  #             mode = "markers", 
  #             text = names, 
  #             hoverinfo = "text",
  #             marker = list(color = get_plat_colors(R, platforms), 
  #                           size = 40))
  toc()
  network
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
  mat <- R$annos[R$clusts[[i]],]
  plots <- c()
  for (n in 1:length(colnames(mat))){
    freqs<- data.frame(table(mat[,n]))
    names(freqs) <- c(colnames(mat)[n], "count")
    plots<- c(plots, ggplot(freqs, aes(x=freqs[,1],y=count,fill =freqs[,1]))+
      geom_bar(stat="identity")+
      coord_flip() +
      xlab(colnames(mat)[n])+
      theme(panel.background = element_blank(), panel.grid.major=element_blank(), legend.position="none")+
      geom_text(aes(x = freqs[,1], y = count, label = count))+
      scale_x_discrete(labels=add_line_format(freqs[,1])))
  }
  plots
}

add_line_format <- function(names){
  unlist(lapply(names, function(i) paste(strwrap(i, 25), collapse="\n")))
}
