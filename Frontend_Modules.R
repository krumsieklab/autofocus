### Frontend Modules ###
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(shinydashboard))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(networkD3))

#### Make Dendrogram ####

#### Label nodes ####
#' get_node_label
#'
#' Gets the label of a node
#' the label is the molecule name if the node is a leaf
#' the label is cluster number, BIC and p-value if the node is internal
#'
#' @param R R struct
#' @param i index of module to be labeled
#'
#' @return label of node i
#'
get_node_label <- function(
    R,
    i
){
  internal_nodes <- dim(R$HCL$merge)[1]

  # Leaf Case
  if (i > (internal_nodes)) return(R$HCL$labels[(i - internal_nodes)] )

  #Internal node case
  else return(paste("<br>Cluster ID:",
                    R$clust_info$ClusterID[[i]] ,
                    "</br><br>Significant Children Density: ",
                    R$clust_info$densities[[i]],
                    "</br><br>Number of Significant Children: ",
                    sum(R$clust_info$densities[(R$clusts[[i]]+dim(R$HCL$merge)[1])])))
}

#' get_node_color
#'
#' Gets the color of a node dependent on user determined threshold
#' the color is the red if module is a significant leaf, black otherwise
#' the color is the red if module is has a higher density than input, black otherwise
#'
#' @param R R struct
#' @param threshold threshold of significant children
#'
#' @return color list
#'
get_node_color <- function(
    R,
    threshold
){
  hc <- R$HCL
  color_list <- rep("black", times = dim(R$clust_info)[1])
  potential_peaks <- identify_peaks(R, dim(hc$merge)[1], threshold, c())
  color_list[potential_peaks]<-'yellow'
  filtered_peaks <- filter_peaks(R, potential_peaks, threshold)
  color_list[filtered_peaks]<-"red"
  sig_kids <- intersect(which(R$clust_info$Size == 1), which(R$clust_info$densities==1))
  color_list[sig_kids]<-"red"
  color_list
}

#' peak_finder_wrapper
#'
#' Wrapper for identify_peaks function
#'
#'
#' @param R R struct
#' @param threshold threshold of significant children
#'
#' @return color list
#'
peak_finder_wrapper<-function (R, threshold){
  hc <- R$HCL
  potential_peaks <- identify_peaks(R, dim(hc$merge)[1], threshold, c())
  filtered_peaks <- filter_peaks(R, potential_peaks, threshold)
  potential_peaks <- potential_peaks[!(potential_peaks %in% filtered_peaks)]
  piggy_backers <- c()
  while(length(potential_peaks)!=0){
    new_potential_peaks <- c()
    for (i in potential_peaks){
      new_potential_peaks <- c(new_potential_peaks,identify_peaks(R, get_disappointment_child(i, R, threshold), threshold, c()))
    }
    new_filtered_peaks <- filter_peaks(R, new_potential_peaks, threshold)
    piggy_backers <- c(piggy_backers, potential_peaks)
    potential_peaks <- new_potential_peaks
    potential_peaks <- potential_peaks[!(potential_peaks %in% new_filtered_peaks)]
    filtered_peaks <- c(filtered_peaks, new_filtered_peaks)
  }

  color_list <- rep("black", times = dim(R$clust_info)[1])
  color_list[piggy_backers]<-'yellow'
  color_list[filtered_peaks]<-"red"
  sig_kids <- intersect(which(R$clust_info$Size == 1), which(R$clust_info$densities==1))
  color_list[sig_kids]<-"red"
  color_list
}

#' identify_peaks
#'
#' Find top peaks at user-defined threshold
#'
#' @param R R struct
#' @param i internal node index to test
#' @param threshold user-defined enrichment threshold
#' @param indices vector collecting peaks
#'
#' @return list of peaks to be considered by user
#'
identify_peaks <- function(
    R,
    i,
    threshold,
    indices
){
  # Leaf Case
  if (i < 0) return(indices)
  else{
    if (R$clust_info$densities[i] > threshold){
      return(c(indices, i))
    }
    else{
      left_peaks <- identify_peaks(R, R$HCL$merge[i,1], threshold, indices)
      right_peaks <- identify_peaks(R, R$HCL$merge[i,2], threshold, indices)
      indices <- c(indices, left_peaks, right_peaks)
    }
  }
  indices
}

#' filter_peaks
#'
#' Filter peaks to get rid of piggy-backers
#'
#' @param R R struct
#' @param peak_list list of potential peaks
#' @param threshold user-defined threshold
#'
#' @return list of peaks to be considered by user
#'
filter_peaks <- function(
    R,
    peak_list,
    threshold
){
  lapply(peak_list, function(i){
    kids <- R$HCL$merge[i,] %>%
      lapply(function(i) if (i < 0) return(abs(i)+dim(R$HCL$merge)[1]) else return(i)) %>% unlist()
    kids_dens <-  R$clust_info$densities[kids]
    if (all(kids_dens>= threshold)){
      return (i)
    }

    else{
      pride <- kids[which(kids_dens == max(kids_dens))]
      questionable_kid <- kids[which(kids_dens == min(kids_dens))]
      if (pride > dim(R$HCL$merge)[1]){
        num_sig <- 1
        length_pride <- 1
      }
      else{
        num_sig <- R$clust_info$densities[sapply(R$clusts[[pride]], function(i) i + dim(R$HCL$merge)[1])] %>% sum()
        length_pride <- R$clusts[[pride]] %>% length()
      }

      length_question <- if (questionable_kid > dim(R$HCL$merge)[1]) 1 else R$clusts[[questionable_kid]] %>% length()

      new_thresh <- num_sig/ (length_pride + length_question)

      if (new_thresh < threshold){
        return(i)
      }
      else{
        return(pride)
      }
    }

  }) %>% unlist()
}

#' get_disappointment_child
#'
#' Determines which child of a node does not contribute
#' to that node's significance density
#'
#' @param i Index of node whose children need to be checked
#' @param R R struct
#' @param threshold threshold of significant children
#'
#' @return The node ID of the child with the lesser density
#'
get_disappointment_child <- function(i, R, threshold){
  kids <- R$HCL$merge[i,] %>%
    lapply(function(i) if (i < 0) return(abs(i)+dim(R$HCL$merge)[1]) else return(i)) %>% unlist()
  kids_dens <-  R$clust_info$densities[kids]
  R$HCL$merge[i,][which(kids_dens == min(kids_dens))]
}

#' get_ancestors
#'
#' Get a list of all internal nodes that are ancestors of input node
#'
#' @param R R struct
#' @param i node id of cluster in question
#' @param parent_list list that is iteratively filled with ancestors
#'
#' @return list of cluster ids of all ancestors of input node
#'
get_ancestors <- function(R, i, parent_list ){
  if (i == dim(R$HCL$merge)[1]){
    return(parent_list)
  }
  else{
    parent <-  which(rowSums(R$HCL$merge==i) == 1)
    parent_list <- c(parent_list, parent)
    get_ancestors(R, parent, parent_list)
  }
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
    R,
    sig_threshold
){
  R$clust_info$label <-mapply(function(i) get_node_label(R, i), 1:nnodes(R$HCL))
  R$clust_info$color <- peak_finder_wrapper(R, sig_threshold)
  colors <- which(R$clust_info$color !="black")
  # Start plotly
  R$clust_info %>% plot_ly(x = ~Coord_X,
                           y = ~Coord_Y, source="A") %>%
    # Add edges between nodes
    add_segments(
      x=R$dend_data$segments$x,
      xend=R$dend_data$segments$xend,
      y=R$dend_data$segments$y,
      yend=R$dend_data$segments$yend)%>%
    # Add node markers
    add_trace(type='scatter',
              mode = "markers",
              marker = list(color=~color)) %>%
    add_trace(x = R$clust_info$Coord_X[colors],
              y = R$clust_info$Coord_Y[colors],
              type='scatter',
              mode = "markers",
              text = R$clust_info$label[colors],
              hovertemplate = paste('<b>%{text}</b>'),
              marker = list(color=R$clust_info$color[colors]))
}


#### Format annotation names ####
#' add_line_format
#'
#' Formats annotation labels adding a line break if they are too long
#'
#' @param names Vector of annotation names
#'
#' @return Vector of annotation names with line split
#'
add_line_format <- function(names){
  unlist(lapply(names, function(i) paste(strwrap(i, 25), collapse="\n")))
}


#' get_drivers
#'
#' Gets the names of the drivers of an mgm
#'
#' @param R R struct
#' @param i cluster index
#'
#' @return names of drivers ofthat cluster's mgm
#'
get_drivers<- function(R, i){
  g <- R$graphs[[i]]
  V(g)$name[V(g)$Driver]
}


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

