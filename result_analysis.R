source(codes.makepath("autofocus/internal_funcs.R"))
library(plotly)

get_sig_nodes <- function(R){
  sig_nodes <- c()
  sig_plats <- c()
  sig_cols <- c()
  for (i in 1:length(R$colors)){
    if (R$colors[i] != "black"){
      sig_nodes <- c(sig_nodes, R$annos$name[unlist(R$clusts[i])])
      sig_plats <- c(sig_plats, R$annos$platform[unlist(R$clusts[i])])
      if (R$colors[i] == "green"){
        sig_cols <- c(sig_cols, R$annos$name[unlist(R$clusts[i])])
      }
    }
  } 
  nodes <- unique(data.frame(sig_nodes, sig_plats))
  nodes$sig_cols <- unlist(lapply(nodes$sig_nodes, function(x) if (x %in% sig_cols) {return('green')} else{return ("NA")}))
  nodes
}


get_sig_clusts <- function(R){
  cluster_no <- c()
  nodes <- c()
  platforms <- c()
  omics <- c()
  sig_cols <- c()
  ps<-c()
  for (i in 1:length(R$colors)){
    if (R$colors[i] != "black"){
      mems <- R$annos$name[unlist(R$clusts[i])]
      ps<-c(ps,rep(R$pvals[i],length(mems)))
      cluster_no <- c(cluster_no, rep(i, length(mems)))
      sig_cols <- c(sig_cols, rep(R$colors[i], length(mems)))
      nodes <- c(nodes, mems)
      platforms <- c(platforms, R$annos$platform[unlist(R$clusts[i])])
      #omics <- c(omics, unlist(lapply(R$annos$platform[unlist(R$clusts[i])], function(x) platform_to_omic(x) )))
    }
  } 
  data.frame(cluster_no, nodes, platforms,  sig_cols,ps)
  
}

platform_to_omic <- function(plat){
  if (plat %in% c("PM","UM","SM","HDF","CM","BM","BRAIN")){
    return("metabolite")
  }
  if (plat == "SOMA"){
    return("protein")
  }
  if (plat %in% c("rawIgG", "rawGP", "IgA")){
    return("glycan")
  }
  if (plat == "LD"){
    return("lipid")
  }
}

sig_nodes <- get_sig_nodes(AD_results)
clusts <- get_sig_clusts(AD_results)

cluster_nos <- unique(clusts$cluster_no)
most_sig_nos <- unique(clusts$cluster_no[which(clusts$sig_cols == "green")])
size <- c()
cross_platform <- c()
cross_omic <- c()
multi_platform_nos <- c()
multi_omic_nos <- c()
for (i in 1:length(cluster_nos)){
  mems <- which(clusts$cluster_no == cluster_nos[i])
  size <- c(size, length(mems))
  cross_platform <- c(cross_platform, length(unique(clusts$platforms[mems])) )
  if (length(unique(clusts$platforms[mems])) > 1){
    multi_platform_nos <- c(multi_platform_nos,cluster_nos[i])
  }
  cross_omic <- c(cross_omic, length(unique(clusts$omics[mems])) )
  if (length(unique(clusts$omics[mems])) > 1){
    multi_omic_nos <- c(multi_omic_nos, cluster_nos[i])
  }
}
plat_freqs <- data.frame(table(sig_nodes$sig_plats))
plat_freqs$most_sig <- table(sig_nodes$sig_plats[which(sig_nodes$sig_cols == "green")])
plat_freqs$tot <- unlist(lapply(plat_freqs$Var1, function(x) sum(rowData(all_platforms)$platform == x)))
pie_fig<- plot_ly(plat_freqs, labels = ~names, values = ~Freq, type = 'pie', textinfo = 'label+value')
xax <- list(
  title = "Platform"
)
yax <- list(
  title = "Number of nodes")
bar_fig <- 
  plot_ly(plat_freqs,
          x = ~names,
          y = ~tot,
          type = 'bar', name = "Number Measured", alpha = 0.5) %>%
  add_trace(plat_freqs,
            x = ~names,
            y = ~Freq,
            type = 'bar', width = 0.6, name = "Number in clusters") %>%
  add_trace(plat_freqs,
            x = ~names,
            y = ~most_sig,
            type = 'bar', width = 0.4, name = "Number in most significant clusters") %>%
  layout(barmode = 'overlay', title = "Platform Representation", xaxis = xax, yaxis = yax)



mod_labels <- c("All Modules", "Tier 1", "Tier 2", "Multi-platform", "Single Platform", "Multi-omic", "Single-omic","Multi-platform", "Single Platform", "Multi-omic", "Single-omic")
mod_colors <- c("green", "orange", 'yellow', "orange", "purple", "red", "blue", "orange","purple","red", "blue")
mod_sources <- c(0,0,1,1,2,2,3,3,7,7)
mod_targets <- c(1,2,3,4,7,8,5,6,9,10)
mod_values <- c(68,188,34,34,135,53,11,23,29,106)
univariate_analysis <- function(
  SE,
  phenotype,
  confounders
){
  
  pvals <- c()
  for (i in 1:nrow(assay(SE))){
    linmod <- lm(colData(SE)[[phenotype]]~assay(SE)[i,] + as.matrix(colData(SE)[,confounders]) + 1)
    pvals <- c(pvals, summary(linmod)$coefficients[2,4])
  }
  rowData(SE)$name[which(pvals < 0.05/length(pvals))]

}


hmdb_overlap <- function(
  SE1, 
  SE2
){
  hmdb1 <- SE1 %>% rowData() %>% as.data.frame() %>% filter(HMDB != "NA" & !is.na(HMDB)) %>% pull(HMDB) %>% remove_bars() %>% unique()
  hmdb2 <- SE2 %>% rowData() %>% as.data.frame() %>% filter(HMDB != "NA" & !is.na(HMDB)) %>% pull(HMDB) %>% remove_bars() %>% unique()
  overlap <- intersect(hmdb1, hmdb2)
  names1 <- lapply(overlap, function(i) {SE1 %>% rowData() %>% as.data.frame() %>% filter(grepl(i, HMDB)) %>% pull(name) %>% paste(sep = "|")}) %>% unlist()
  names2 <- lapply(overlap, function(i) {SE2 %>% rowData() %>% as.data.frame() %>% filter(grepl(i, HMDB)) %>% pull(name) %>% paste(sep = "|")}) %>% unlist()
  cbind(overlap, names1, names2)
}

remove_bars <- function(hmdb_list){
  unlist(lapply(hmdb_list, function(i) unlist(strsplit(i, split='|', fixed=TRUE))))
}



sankey_diag <- function(
  labels,
  colors,
  values,
  sources,
  targets
){

  
  fig <- plot_ly(
    type = "sankey",
    orientation = "h",
    
    node = list(
      label = labels,#c("Significant hits", "Leaves", "Modules", "Found in Univariate", "New Hits", "Tier 1", "Tier 2"),
      color = colors,#c("green", "blue", "yellow", "blue", "purple", "orange", "yellow"),
      pad = 15,
      thickness = 20,
      line = list(
        color = "black",
        width = 0.5
      )
    ),
    
    link = list(
      source = sources,#c(0,1,1,0,2,2),
      target = targets,#c(1,4,3,2,6,5),
      value =  values#c(90,10,80,256,188,68)
    )
  )
  fig <- fig %>% layout(
    title = "Significant Hits Composition",
    font = list(
      size = 16
    )
  )
  
  fig
}


vline <- function(x = 0, color = "red") {
  list(
    type = "line", 
    y0 = 0, 
    y1 = 2000, 
    yref = "paper",
    x0 = x, 
    x1 = x, 
    line = list(color = color)
  )
}
cor_mat <- cor(t(assay(all_platforms)))
hmdbs <- rowData(all_platforms)$HMDB %>% unique() %>% strsplit("\\|") %>% unlist()
hmdbs <- hmdbs[!(is.na(hmdbs))]
hmdbs  <- hmdbs[which(hmdbs!="NA")]
cors <- c()
for (i in 1:length(hmdbs)){
  sames <- which(grepl(rowData(all_platforms)$HMDB, hmdbs[i]))
  if (length(sames)>1){
    same_cors <- cor_mat[sames, sames]
    same_cors <- same_cors[upper.tri(same_cors)]
    cors <- c(cors, same_cors)
  }
}
  fig <- plot_ly(alpha = 0.6) %>% 
    add_histogram(x = cors) %>% 
    layout(title = "Histogram of same molecule correlations")

  fig

  
y_hat_y <- function(
  R,
  confounders
){
  cors <- c()
  cols <- c()
   for(i in 1:length(R$colors)) {
    if (i < 5135){
      mems<- R$annos$name[unlist(R$clusts[i])]
      names <- rowData(all_platforms)$FeatureId[which(rowData(all_platforms)$name %in% mems)]
      full_data <- cbind(t(assay(all_platforms)),as.matrix(R$samples[,confounders]))
      model <- lm(colData(all_platforms)$DIAB ~ full_data[,names] + full_data[,confounders]+1)
      coeffs <- summary(model)$coefficients[,"Estimate"]
      partial_data <- cbind(rep(1, dim(all_platforms)[2]),full_data[,names], full_data[,confounders])
      if (length(coeffs) == dim(partial_data)[2] & length(mems) < 51){
        y_hat <- partial_data %*% coeffs
        cors <- c(cors, cor(y_hat, colData(all_platforms)$DIAB))
        cols <- c(cols, R$colors[i])
      }
      }
  }
  data.frame(cors, cols)
}


better_net <- function(
  R, 
  i
){
  members <- R$clusts[[3814]]
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
  rownames(adj_mat) <- names
  colnames(adj_mat)<- names
  diag(adj_mat)<-0
  G<- graph_from_adjacency_matrix(adj_mat)
  test.visn <- toVisNetworkData(G)
  test.visn$edges$value <- G$edges$weight
  test.visn$edges$color = "black"
  test.visn$nodes$color =  c("red", "blue", "red", "red", "red" )
  test.visn$nodes$shape = c("triangle", "square","triangle","triangle","triangle")
  test.visn$nodes$label = NULL
  visNetwork(test.visn$nodes, test.visn$edges)
  
}





