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
  if("densities" %in% names(R$clust_info)){
    return(paste("<br>Cluster ID:",
                 R$clust_info$ClusterID[[i]] ,
                 "</br><br>Significant Children Density: ",
                 R$clust_info$densities[[i]],
                 "</br><br>Number of Significant Children: ",
                 sum(R$clust_info$densities[(R$clusts[[i]]+dim(R$HCL$merge)[1])])))
  }else{
    return(paste("<br>Cluster ID:",
                 R$clust_info$ClusterID[[i]] ,
                 "</br><br>Significant Children Density - Phenotype 1: ",
                 R$clust_info$pheno1_densities[[i]],
                 "</br><br>Number of Significant Children - Phenotype 1: ",
                 sum(R$clust_info$pheno1_densities[(R$clusts[[i]]+dim(R$HCL$merge)[1])]),
                 "</br><br>Significant Children Density - Phenotype 2: ",
                 R$clust_info$pheno2_densities[[i]],
                 "</br><br>Number of Significant Children - Phenotype 2: ",
                 sum(R$clust_info$pheno2_densities[(R$clusts[[i]]+dim(R$HCL$merge)[1])])
    ))
  }

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
  if("densities" %in% names(R$clust_info)){
    potential_peaks <- identify_peaks(R, dim(hc$merge)[1], threshold, c())
    filtered_peaks <- filter_peaks(R, potential_peaks, threshold)
    potential_peaks <- potential_peaks[!(potential_peaks %in% filtered_peaks)]
    # piggy_backers <- c()
    # while(length(potential_peaks)!=0){
    #   new_potential_peaks <- c()
    #   for (i in potential_peaks){
    #     new_potential_peaks <- c(new_potential_peaks,identify_peaks(R, get_disappointment_child(i, R, threshold), threshold, c()))
    #   }
    #   new_filtered_peaks <- filter_peaks(R, new_potential_peaks, threshold)
    #   piggy_backers <- c(piggy_backers, potential_peaks)
    #   potential_peaks <- new_potential_peaks
    #   potential_peaks <- potential_peaks[!(potential_peaks %in% new_filtered_peaks)]
    #   filtered_peaks <- c(filtered_peaks, new_filtered_peaks)
    # }

    color_list <- rep("#666666", times = dim(R$clust_info)[1])
    #color_list[piggy_backers]<-'yellow'
    color_list[filtered_peaks]<-"#E7298A"
    sig_kids <- intersect(which(R$clust_info$Size == 1), which(R$clust_info$densities==1))
    color_list[sig_kids]<-"#E7298A"
    color_list

  }else if("pheno1_densities" %in% names(R$clust_info) && "pheno2_densities" %in% names(R$clust_info)){
    # get peaks and piggybackers for phenotype 1
    potential_peaks <- identify_peaks(R, dim(hc$merge)[1], threshold, c(), "pheno1_densities")
    filtered_peaks <- filter_peaks(R, potential_peaks, threshold, "pheno1_densities")
    potential_peaks <- potential_peaks[!(potential_peaks %in% filtered_peaks)]
    piggy_backers <- c()
    while(length(potential_peaks)!=0){
      new_potential_peaks <- c()
      for (i in potential_peaks){
        new_potential_peaks <- c(new_potential_peaks,identify_peaks(R, get_disappointment_child(i, R, threshold, "pheno1_densities"), threshold, c(), "pheno1_densities"))
      }
      new_filtered_peaks <- filter_peaks(R, new_potential_peaks, threshold)
      piggy_backers <- c(piggy_backers, potential_peaks)
      potential_peaks <- new_potential_peaks
      potential_peaks <- potential_peaks[!(potential_peaks %in% new_filtered_peaks)]
      filtered_peaks <- c(filtered_peaks, new_filtered_peaks)
    }
    pheno1_peaks <- filtered_peaks
    # pheno1_piggy <- piggy_backers

    # get peaks and piggybackers for phenotype 2
    potential_peaks <- identify_peaks(R, dim(hc$merge)[1], threshold, c(), "pheno2_densities")
    filtered_peaks <- filter_peaks(R, potential_peaks, threshold, "pheno2_densities")
    potential_peaks <- potential_peaks[!(potential_peaks %in% filtered_peaks)]
    piggy_backers <- c()
    while(length(potential_peaks)!=0){
      new_potential_peaks <- c()
      for (i in potential_peaks){
        new_potential_peaks <- c(new_potential_peaks,identify_peaks(R, get_disappointment_child(i, R, threshold, "pheno2_densities"), threshold, c(), "pheno2_densities"))
      }
      new_filtered_peaks <- filter_peaks(R, new_potential_peaks, threshold)
      piggy_backers <- c(piggy_backers, potential_peaks)
      potential_peaks <- new_potential_peaks
      potential_peaks <- potential_peaks[!(potential_peaks %in% new_filtered_peaks)]
      filtered_peaks <- c(filtered_peaks, new_filtered_peaks)
    }
    pheno2_peaks <- filtered_peaks
    # pheno2_piggy <- piggy_backers

    # find unique and overlapping colors for the two phenotypes
    common_peaks <- intersect(pheno1_peaks, pheno2_peaks)
    # common_piggy <- intersect(pheno1_piggy, pheno2_piggy)

    color_list <- rep("#666666", times = dim(R$clust_info)[1])
    color_list[pheno1_peaks]<-'#E7298A'
    color_list[pheno2_peaks]<-"#A6D854"
    color_list[common_peaks] <- "#FC8D62"
    # color_list[pheno1_piggy] <- "orange"
    # color_list[pheno2_piggy] <- "green"
    # color_list[common_piggy] <- "yellow"

    # get color of significant leaves
    pheno1_sig_kids <- intersect(which(R$clust_info$Size == 1), which(R$clust_info$pheno1_densities==1))
    pheno2_sig_kids <- intersect(which(R$clust_info$Size == 1), which(R$clust_info$pheno2_densities==1))
    common_sig_kids <- intersect(pheno1_sig_kids, pheno2_sig_kids)
    color_list[pheno1_sig_kids]<-"#E7298A"
    color_list[pheno2_sig_kids] <- "#A6D854"
    color_list[common_sig_kids] <- "#FC8D62"
    color_list

  }else{
    stop("R$clust_info contains no valid column names for densities")
  }
}

#' Find top peaks at user-defined threshold
#'
#' ADD DESCRIPTION.
#'
#' @param R R struct
#' @param i internal node index to test
#' @param threshold user-defined enrichment threshold
#' @param indices Vector collecting peaks
#' @param dense_col Name of the column in the clust_info data frame containing densities for the
#'    phenotype of interest.
#'
#' @return list of peaks to be considered by user
#'
#' @noRd
identify_peaks <- function(R, i, threshold, indices, dense_col="densities"){
  # Leaf Case
  if (i < 0) return(indices)
  else{
    if (R$clust_info[[dense_col]][i] >= threshold){
      return(c(indices, i))
    }
    else{
      left_peaks <- identify_peaks(R, R$HCL$merge[i,1], threshold, indices, dense_col)
      right_peaks <- identify_peaks(R, R$HCL$merge[i,2], threshold, indices, dense_col)
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
#' @param dense_col Name of column in clust_info data frame that contains the densities for the
#'    phenotype of interest.
#'
#' @return list of peaks to be considered by user
#'
#' @noRd
filter_peaks <- function(R, peak_list, threshold, dense_col = "densities"){

  lapply(peak_list, function(i){
    kids <- R$HCL$merge[i,] %>%
      lapply(function(i) if (i < 0) return(abs(i)+dim(R$HCL$merge)[1]) else return(i)) %>% unlist()
    kids_dens <-  R$clust_info[[dense_col]][kids]

    if (all(kids_dens>= threshold)){
      return (i)
    }else{
      pride <- kids[which(kids_dens == max(kids_dens))]
      questionable_kid <- kids[which(kids_dens == min(kids_dens))]
      if (pride > dim(R$HCL$merge)[1]){
        num_sig <- 1
        length_pride <- 1
      }
      else{
        num_sig <- R$clust_info[[dense_col]][pride] * R$clust_info$Size[pride]
        length_pride <- R$clusts[[pride]] %>% length()
      }

      num_sig_questionable <- R$clust_info[[dense_col]][questionable_kid] * R$clust_info$Size[questionable_kid]
      if (num_sig_questionable ==0){
        return(pride)
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
get_disappointment_child <- function(i, R, threshold, dense_col="densities"){
  kids <- R$HCL$merge[i,] %>%
    lapply(function(i) if (i < 0) return(abs(i)+dim(R$HCL$merge)[1]) else return(i)) %>% unlist()
  kids_dens <-  R$clust_info[[dense_col]][kids]
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
  R$clust_info$label <-mapply(function(i) get_node_label(R, i), 1:dendextend::nnodes(R$HCL))
  R$clust_info$color <- peak_finder_wrapper(R, sig_threshold)
  colors <- which(R$clust_info$color !="#666666")
  # Start plotly
  R$clust_info %>% plotly::plot_ly(x = ~Coord_X,
                                   y = ~Coord_Y, source="A") %>%
    # Add edges between nodes
    plotly::add_segments(
      x=R$dend_data$segments$x,
      xend=R$dend_data$segments$xend,
      y=R$dend_data$segments$y,
      yend=R$dend_data$segments$yend)%>%
    # Add node markers
    plotly::add_trace(type='scatter',
                      mode = "markers",
                      text = R$clust_info$label,
                      hovertemplate = paste('<b>%{text}</b>'),
                      marker = list(color=~color)) %>%
    plotly::add_trace(x = R$clust_info$Coord_X[colors],
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
  if(length(g)==1) return(character(0))
  igraph::V(g)$name[igraph::V(g)$Driver]
}

#' Get pie plot colors
#'
#' Assign colors to each unique value of an annotation.
#'
#' @param anno_vec Vector of unique annotation values.
#'
#' @return List of unique annotation values with colors as names.
#'
#' @importFrom RColorBrewer brewer.pal
#'
#' @noRd
get_pie_color <- function(anno_vec){

  palette <- unique(c(brewer.pal(n = 12,name="Set3"),
                      brewer.pal(n = 9,name="Set1"),
                      brewer.pal(n = 8,name="Set2"),
                      brewer.pal(n = 9,name="Pastel1"),
                      brewer.pal(n = 8,name="Pastel2"),
                      brewer.pal(n = 12,name="Paired"),
                      brewer.pal(n = 8,name="Dark2")))
  palette <- sample(palette)

  anno_vals <- unique(anno_vec)
  if(length(anno_vals) > length(palette)) stop("Too many unique values for annotation plot.")

  anno_vals <- anno_vals[anno_vals != "No Data"]
  pal <- as.list(c("gray", palette[1:length(anno_vals)]))
  anno_vals <- c("No Data", anno_vals)
  names(pal) <- anno_vals

  as.list(pal)
}

#' Template pie plot function
#'
#' This is a template function called to dynamically define pie plot columns in the reactable
#' table. The object 'anno' (surrounded by {{ }}) will be replaced with the annotation column name
#' of interest using the function tmpl(). This function defintion is then provided to the function
#' colDef() via the cell argument (see tmpl and colDef documentation for details).
#'
#' @param value The cell info object. See reactable::colDef for details.
#'
#' @return Function definition for the custum cell renderer. See reactable::colDef for details.
#'
#' @noRd
pie_plot_tmpl_fun <- function(value){

  pie_df <- anno_df %>%
    filter(cluster == value) %>%
    group_by(cluster, {{anno}}) %>%
    summarise(n=n()) %>%
    mutate(freq=n/sum(n)) %>%
    ungroup()
  pal <- pal_list %>% extract2("{{anno}}")
  pie_df$colors <- unname(pal[match(unlist(dplyr::select(pie_df, {{anno}})),names(pal))])
  #color_vec <- pal[match(unlist(dplyr::select(pie_df, {{anno}})),names(pal))]
  p <- plot_ly(pie_df,
               labels = ~{{anno}},
               values = ~freq,
               type="pie",
               hoverinfo='label',
               marker=list(colors = ~colors),
               rotation=300,
               textfont = list(size = 20,color="black")) %>%
    #marker=list(colors=names(color_vec[which(color_vec %in% unlist(dplyr::select(pie_df, SUPER_PATHWAY)))]))) %>%
    layout()

  return(p)

}
