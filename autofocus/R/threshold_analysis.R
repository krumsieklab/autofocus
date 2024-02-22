#' Perform threshold analysis 
#'
#' Outputs a plot across thresholds from 0.01 to 0.99 at 0.01 increments with the following information:
#' \itemize{
#'    \item{Cluster Height Range: The height range of returned clusters at that threshold, normalized to the largest range across all thresholds}
#'    \item{Cluster Size Range: The size range of returned clusters at that threshold, normalized to the largest range across all thresholds}
#'    \item{Cluster Count: The number of returned clusters at that threshold, normalized to the highest number of clusters across all thresholds}
#'    \item{Non-Significant Node Proportion: The proportion of nodes in returned clusters that are not significantly associated with the phenotype of interest, normalized to the largest proportion across all thresholds}
#'    \item{Singleton Node Count: The number of nodes significantly associated to phenotype that are not found in a cluster, normalized to the largest number of singletons across all thresholds}
#'    \item{Noise Minimizing Score: A metric that combines the non-significant node proportion and singleton node count, to show at which threshold they are minimized}
#' }
#'
#' @param R R_struct object output by initialize_R
#' @param phenotype The phenotype for which to perform the threshold analysis (clusters are phenotype dependent)
#' @param phenotype_name Formal name of phenotype for plotting title
#' 
#' @return None
#'
#'
#' @author AS
#'
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#' @import RColorBrewer
#'
#' @export
threshold_analysis <- function(R, phenotype, phenotype_name = NA){
  
  if(!(phenotype %in% R$phenotypes)){
    stop("Phenotype not found in R struct")
  }
  
  if (length(R$phenotypes) >1){
    if (R$phenotypes[1] == phenotype){
      R$clust_info$densities <- R$clust_info$pheno1_densities
    }
    else{
      R$clust_info$densities <- R$clust_info$pheno2_densities
    }
  }
  
  heights_to_test <- seq(0.01,0.99,0.01)
  
  threshold_range_results <- lapply(heights_to_test, function (i) {
    height_range <- get_cluster_info(R, i , metric = "height_range")
    size_range <- get_cluster_info(R, i , metric = "clust_size_range")
    prop_non_sig_to_sig <- get_cluster_info(R, i)
    num_clusts <- get_cluster_info(R, i, metric="num_clusts")
    singletons <- get_cluster_info(R, i, metric = "singletons")
    data.frame(threshold=i, 
               height_range, 
               size_range, 
               prop_non_sig_to_sig, 
               num_clusts,
               singletons)
  }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(height_range_norm = height_range/max(height_range)) %>% 
    dplyr::mutate(prop_non_sig_norm = prop_non_sig_to_sig/max(prop_non_sig_to_sig)) %>% 
    dplyr::mutate(singleton_norm = singletons/max(singletons)) %>%
    dplyr::mutate(num_clusts_norm = num_clusts/max(num_clusts)) %>% 
    dplyr::mutate(clust_size_range_norm = size_range/max(size_range)) %>% 
    dplyr::mutate(score = prop_non_sig_norm + singleton_norm) %>% 
    dplyr::mutate(score_norm  = min(score)/score) 
  
  
  disc_palette <- c(brewer.pal(n=8, name= "Dark2"), brewer.pal(n=8, name= "Set2"))
  
  if(is.na(phenotype_name)) phenotype_title=phenotype else phenotype_title=phenotype_name
  
  threshold_range_results %>% 
    dplyr::select(threshold, 
                  height_range_norm,
                  clust_size_range_norm,
                  num_clusts_norm,
                  prop_non_sig_norm, 
                  singleton_norm, 
                  score_norm) %>% 
    reshape2::melt(id="threshold") %>% 
    ggplot(aes(x=threshold, y=value, color=variable))+
    geom_line()+ 
    ggtitle(paste0("Threshold Analysis Results for ", phenotype_title))+
    scale_color_manual(values=disc_palette[1:6], 
                       name= "Cluster Metric",
                       labels = c("Cluster Height Range", 
                                  "Cluster Size Range",
                                  "Cluster Count",
                                  "Non-Significant Node Proportion",
                                  "Singleton Node Count",
                                  "Noise Minimizing Score")) +
    theme_minimal()+
    geom_label(data = threshold_range_results %>%
                       select(threshold, score_norm) %>%
                       filter(score_norm == max(score_norm)) %>%
                       melt(id="threshold") %>%
                       .[1,],
               aes(x = threshold,
                        y = value,
                        label = paste0("Noise minimized at threshold: ", threshold)),
               show.legend=F,
               vjust = -0.5) +
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10,face="bold"))+
  geom_point() +
   xlab("Threshold")+
   ylab("Normalized metric value")+
   guides(color = guide_legend(override.aes = list(size = 5)))
  
}


