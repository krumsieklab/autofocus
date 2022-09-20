#' Get node label
#'
#' Gets the label of a node. The label is the molecule name if the node is a leaf. The label is
#' the cluster number, BIC and p-value if the node is internal.
#'
#' @param R R struct object.
#' @param i The index of the module to be labeled.
#'
#' @return Label of node i.
#'
#' @noRd
get_node_label <- function(R, i){
  internal_nodes <- dim(R$HCL$merge)[1]

  # Leaf Case
  if (i > (internal_nodes)) return(R$HCL$labels[(i - internal_nodes)] )

  #Internal node case
  return(paste("<br>Cluster ID:",
               R$clust_info$ClusterID[[i]] ,
               "</br><br>Significant Children Density: ",
               R$clust_info$densities[[i]],
               "</br><br>Number of Significant Children: ",
               sum(R$clust_info$densities[(R$clusts[[i]]+dim(R$HCL$merge)[1])])))
}

#' Wrapper for identify_peaks function
#'
#' ADD DESCRIPTION
#'
#' @param R R struct object
#' @param threshold Threshold of significant children.
#'
#' @return color list
#'
#' @noRd
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

#' Find top peaks at user-defined threshold
#'
#' ADD DESCRIPTION.
#'
#' @param R R struct
#' @param i internal node index to test
#' @param threshold user-defined enrichment threshold
#' @param indices Vector collecting peaks
#'
#' @return list of peaks to be considered by user
#'
#' @noRd
identify_peaks <- function(R, i, threshold, indices){
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

#' Filter peaks
#'
#' Filter peaks to get rid of piggy-backers.
#'
#' @param R R struct
#' @param peak_list list of potential peaks
#' @param threshold user-defined threshold
#'
#' @return list of peaks to be considered by user
#'
#' @noRd
filter_peaks <- function(R, peak_list, threshold){
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

#' Get 'disappointment' child
#'
#' Determines which child of a node does not contribute to that node's significance density.
#'
#' @param i Index of node whose children need to be checked
#' @param R R struct
#' @param threshold threshold of significant children
#'
#' @return The node ID of the child with the lesser density
#'
#' @noRd
get_disappointment_child <- function(i, R, threshold){
  kids <- R$HCL$merge[i,] %>%
    lapply(function(i) if (i < 0) return(abs(i)+dim(R$HCL$merge)[1]) else return(i)) %>% unlist()
  kids_dens <-  R$clust_info$densities[kids]
  R$HCL$merge[i,][which(kids_dens == min(kids_dens))]
}

#' Get a list of all internal nodes that are ancestors of input node
#'
#' ADD DESCRIPTION
#'
#' @param R R struct
#' @param i node id of cluster in question
#' @param parent_list list that is iteratively filled with ancestors
#'
#' @return list of cluster ids of all ancestors of input node
#'
#' @noRd
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

#' Create dendrogram view of data
#'
#' ADD DESCRIPTION
#'
#' @param R R struct
#' @param order_coords coordinates of points in plot
#'
#' @return dend_network which is a plotly network view of the dendrogram
#'
#' @noRd
plot_dend<-function(R, sig_threshold){
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

#' Get drivers
#'
#' Gets the names of the drivers of an mgm
#'
#' @param R R struct
#' @param i cluster index
#'
#' @return names of drivers ofthat cluster's mgm
#'
#' @noRd
get_drivers<- function(R, i){
  g <- R$graphs[[i]]
  V(g)$name[V(g)$Driver]
}
