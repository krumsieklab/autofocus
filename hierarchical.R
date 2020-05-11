#### init ----
source(codes.makepath("autofocus/internal_funcs.R"))


get_node_color <- function(R, i, signif){
  hc <- R$HCL
  internal_nodes <- dim(R$HCL$merge)[1]
  if (i %in% signif[!is.na(signif)]){
    if (i > (internal_nodes)) return('blue')
    else {
      
      children <- hc$merge[i,]
      child1 <- if (children[1]<0) (internal_nodes + abs(children[1])) else children[1]
      child2 <- if (children[2]<0) (internal_nodes + abs(children[2])) else children[2]
      
      ## Children aren't significant, combination is
      if (R$pvals[child1] > 0.05 & R$pvals[child2] > 0.05) return("purple")
      
      ## Children are significant
      if(R$pvals[child1] < 0.05 & R$pvals[child2] < 0.05){
        
        if (R$BIC[i] < R$BIC[child1] & R$BIC[i]<R$BIC[child2]) return("purple") else return("orange")
      }
      
      ## One child is significant
      else {
        sig_child <- if (R$pvals[child1] < R$pvals[child2]) child1 else child2
        
        if(R$BIC[i] < R$BIC[sig_child]) return("purple") else return("orange")
      }
    }
  }
  else{
    return("black")
  }
}


module.pval <- function(i, full_data, R, phenotype, confounders, return_BIC = FALSE, max_size){
  me <- R$HCL$merge
  if (i <= dim(me)[1]){
    members <- R$clusts[i][[1]]
    if (length(members)>max_size){
      return(NA)
    }
  }
  else{
    members <- (i - dim(me)[1])
  }
  
  if (length(confounders)==0){
    # test model
    model.test <- lm(phenotype ~ full_data[,members] + 1) 
    # base model
    model.base <- lm(phenotype ~ 1)
  }
  else{
    # test model
    model.test <- lm(phenotype ~ full_data[,members] + full_data[,confounders]+1) 
    # base model
    model.base <- lm(phenotype ~ full_data[,confounders]+1) 
  }
  # p-val of model differences
  if (return_BIC){
    return(BIC(model.test))
  }
  anova(model.test, model.base)$`Pr(>F)`[2]
}

find_sig_clusts <- function(R, outcome, confounders, cores = 1){
  hc <- R$HCL
  
  PHEN <- R$samples[[outcome]]
  X <- R$data
  #### scoring function definition ----
  full_data <- cbind(X, PHEN, as.matrix(R$samples[,confounders]))
  
  #### find all valid clusters and score them ----
  ###CUT 1###
  
  ## calculate p-values
  internal_nodes <- dim(hc$merge)[1]
  leaves <- nleaves(hc)
  all_nodes <- nnodes(hc)
  print("Calculating pvals")
  allpvals <- sapply(1:all_nodes, module.pval, full_data, R, PHEN, confounders, max_size = 50); 
  print("Calculating BIC")
  R$BIC <- sapply(1:all_nodes, module.pval, full_data, R, PHEN, confounders, max_size = 50, TRUE);
  
  inds <- which(!is.na(allpvals))
  
  #### WY p-values ----
  ## -> randomize outcome, run all tests, record smallest p-value of each iterations
  
  # server timing, with old version, not all nodes (singltons missing, real time will be about double)
  
  # - nrand=100000, 40 cores, 135184.517 sec elapsed (37.5h)
  
  # parameters
  nrand <- 1000
  mc.cores <- cores
  # safe parallelization
  
  # run
  
  if (mc.cores>1) doParallel::registerDoParallel(cores=mc.cores)
  wy.null <- foreach(r = 1:nrand) %dopar% {  min(sapply(inds, module.pval, full_data, R, sample(PHEN), confounders, max_size = 50)) } %>% unlist()
  
  # NAs can occur (very rarely), cut them out
  cut <- which(is.na(wy.null))
  if (length(cut)>0) {
    print(sprintf("NAs found in random samples: %s", paste0(cut, collapse = ',')))
    wy.null <- wy.null[-cut]
    nrand <- length(wy.null)
    # plot null distribution
    hist(log10(wy.null))
    plot(log10(wy.null))
  }
  
  ## empirical WY adjusted p-value generation
  # -> count how many of the WY null values are actually lower than the observed p-value
  R$pvals <- sapply(allpvals, function(pval){sum(wy.null<pval)}) / nrand
  ###CUT 2###
  
  # determine significant nodes to be colored
  print("Getting node colors")
  signif <- which(R$pvals<0.05)
  
  R$colors <- mapply(function(i) get_node_color(R, i, signif), 1:all_nodes)
  ###CUT 3###
  R<- R %>% build_sig_clusts() %>% make_full_dend() %>% make_bar_plots() %>% make_windows()
  to_remove <- c("data", "C", "dist", "clusts", "annos")
  R[!(names(R) %in% to_remove)]  
}


# 
# 
# 
# 
# ###CUT 1###
# ## Mustafa tree walker function (including old ones)
# fJan <- function(hc, min_size, max_size){
#   
#   get_uc_vc <-function(hc, th){
#     # # of samples = # of possible clustering including 1 cluster
#     k = length(hc$order)
#     
#     # cutting the tree every possible level possible given th
#     omat <- cutree(hc, k = seq(length(hc$order)-th-1))[hc$order,]
#     omat = omat[,!apply(omat,2, function(x) all(table(x) < th) )]
#     
#     # giving a hashed id to every split at every level as hashed.id = i1_i2
#     chah <- 
#       apply(omat,2, function(x)
#         sapply(unique(x), function(i){
#           inds = which(x==i)
#           paste(min(inds), max(inds), sep = "_") 
#         })
#       ) 
#     
#     # get unique splits on the tree with level
#     csm = sapply(2:length(chah), function(i) c(setdiff(chah[[i]],chah[[i-1]]),i) )
#     rm(chah)
#     # set cluster ids
#     mmat = t(matrix(as.numeric( factor(c(csm[1:2,]), levels = c(csm[1:2,]))), ncol = ncol(csm)))+1
#     # locations of clusters
#     sz = sapply(apply(csm[1:2,],2,strsplit,split="_"), function(x) as.numeric(unlist(x)))
#     l = as.numeric(csm[3, ])
#     mmat = cbind(t(sz), l = l, h = rev(hc$height)[l-1], mmat)
#     # sizes of clusters
#     sz = apply(sz,2,diff)[-2,]+1
#     rm(csm,l)
#     # to keep, at least 1 valid clusters at level h  
#     mmat = mmat[!(colSums(sz<th)==2), ]
#     sz = sz[,!(colSums(sz<th)==2)]
#     colnames(mmat) = c("c1_1", "c1_2", "c2_1", "c2_2", "level", "h", "cid1", "cid2")
#     # unvalid cluster of the pair
#     mmat[sz[1,]<th,"cid1"] = NA
#     mmat[sz[2,]<th,"cid2"] = NA
#     rownames(sz) = c("n_c1","n_c2")
#     # set sized for non-hiertest
#     sz[1,is.na(mmat[,"cid1"])] = k - sz[2,is.na(mmat[,"cid1"])]
#     sz[2,is.na(mmat[,"cid2"])] = k - sz[1,is.na(mmat[,"cid2"])]
#     
#     list( clinds = data.frame(cbind(mmat,t(sz))), minsize = th )
#   }
#   sm = get_uc_vc(hc, min_size)$clinds
#   
#   c1 = sm[,1:2]
#   colnames(c1)=c("i","j")
#   l1 = c1$j-c1$i
#   c1 = c1[l1>=min_size & l1<max_size, ]
#   
#   c2 = sm[,3:4]
#   colnames(c2)=c("i","j")
#   l2 = c2$j-c2$i
#   c2 = c2[l2>=min_size & l2<max_size, ]
#   
#   structure(rbind(c1,c2), min_size=min_size, max_size=max_size)
# }
# # function that solves the problem recursively
# fJan2 <- function(hc, min_size=0, max_size=Inf){
#   
#   fwrap <- function(x,me, min_size, max_size){
#     
#     ff_<-function(x,cl=c()){
#       if(sum(x<0)==2) return(c(cl,-x))
#       if(x[1]<0 & x[2]>0) return(ff_(me[x[2],],c(cl,-x[1])))
#       if(x[2]<0 & x[1]>0) return(ff_(me[x[1],],c(cl,-x[2])))
#       return(c(cl,ff_(me[x[1],]),ff_(me[x[2],])))
#     }
#     
#     cl = ff_(x)
#     k = length(cl)
#     if(min_size>k || k>max_size) return(NULL)
#     cl
#   }
#   
#   lst <- apply(hc$merge, 1, fwrap, me=hc$merge, min_size=min_size, max_size=max_size)
#   # invert
#   lst <- lst[length(lst):1]
#   # add singletons
#   lst <- c(lst, hc$order)
# }
# # dendextend naming compatible tree walker
# fJan3 <- function(hc, min_size=0, max_size=Inf){
#   
#   fwrap <- function(x,me, min_size, max_size){
#     
#     ff_<-function(x,cl=c()){
#       if(sum(x<0)==2) return(c(cl,-x))
#       if(x[1]<0 & x[2]>0) return(ff_(me[x[2],],c(cl,-x[1])))
#       if(x[2]<0 & x[1]>0) return(ff_(me[x[1],],c(cl,-x[2])))
#       return(c(cl,ff_(me[x[1],]),ff_(me[x[2],])))
#     }
#     
#     cl = ff_(x)
#     k = length(cl)
#     if(min_size>k || k>max_size) return(NA)
#     x
#   }
#   
#   rel <- apply(hc$merge, 1, fwrap, me=hc$merge, min_size=min_size, max_size=max_size)
#   
#   # depth first search for naming
#   fDFS <- function(hc){
#     me = hc$merge
#     tx = me*0
#     dfs <- function(i, tb, k =1){
#       if(sum(me[i,]<0)==2){
#         tb[i, ] = c(k+1,k+2)
#         return(tb)
#       }
#       if(me[i,1]<0){
#         tb[i,] = c(k+1,k+2)
#         return(dfs(me[i,2],tb,k+2))
#       }
#       if(me[i,2]<0){
#         tb[i,1] = k+1
#         tb = dfs(me[i,1],tb,k+1)
#         tb[i,2] = max(tb)+1
#         return(tb)
#       }
#       tb[i,1]= k+1
#       tb = dfs(me[i,1],tb, k+1)
#       tb[i,2] = max(tb)+1
#       tb = dfs(me[i,2],tb, max(tb))
#       return(tb)
#     }
#     dfs(nrow(tx),tx)
#   }
#   nmm = fDFS(hc) 
#   me = hc$merge
#   singletons = abs(me[me < 0])
#   names_singletons = nmm[me < 0]
#   if(min_size>1) singletons = rep(list(NA), length(singletons))
#   
#   c(rel[c(me[me > 0], nrow(me))], singletons)[order(c(c(nmm[me > 0],1), names_singletons))]
# }
# 
# # old code for old function versions from Mustafa
# # tic("fJan"); sm <- fJan(hc, min_size = 1, max_size = 50); toc()
# # tic("map back"); inds<-apply(sm,1,function(x) (seq(hc$order)[hc$order])[x[1]:x[2]]); toc()
# # tic("fJan2"); sm2 <- fJan2(hc, min_size = 1, max_size = 50); toc()
# # tic("all lm models on sm1"); allpvals <- sapply(inds, module.pval, full_data[,'T2D']); toc()
# 
# # ## find clusters
# tic("fJan3"); sm3 <- fJan3(hc, min_size = 1, max_size = 50); toc()
# # convert it to a list version with no NAs, and node IDs as names
# nona <- !sapply(sm3, function(x){all(is.na(x))})
# inds <- sm3[nona]
# names(inds) <- nona %>% which() %>% as.character()
# 
# 
# ###CUT 2###
# plot(R$pvals)
# plot(sort(R$pvals))
# 
# ## fit gamma distribution to -log10 pvals
# wy.null.log10 <- -log10(wy.null)
# library(fitdistrplus)
# fit.gamma <- fitdist(wy.null.log10, distr = "gamma", method = "mle")
# summary(fit.gamma); plot(fit.gamma)
# # fit.lnorm <- fitdist(wy.null.log10, distr = "lnorm", method = "mle")
# # summary(fit.lnorm); plot(fit.lnorm)
# 
# 
# ## derive estimated p-values from theoretical distribution
# R$pvals.est <- pgamma(-log10(allpvals), shape=fit.gamma$estimate['shape'], rate = fit.gamma$estimate['rate'], lower.tail = F)
# # compare, but only where empirical adjusted pval is not 0 or 1
# use <- R$pvals>0 & R$pvals<1
# plot(log10(R$pvals[use]), log10(R$pvals.est[use]))
# abline(0,1,col='red')
# 
# 
# ###CUT 3 ###
# # initialize dendrogram visualization
# dendro <- hc %>% 
#   as.dendrogram() %>% 
#   dendextend::set("labels", NA) %>%
#   dendextend::set("nodes_pch", 19) %>%
#   dendextend::set("nodes_col", 4)
# # show all significant clusters
# dendro.sign <- dendro %>% dendextend::set("nodes_cex", marksignif) 
# dendro.sign %>% plot(main='all significant clusters')
# # highlight some examples (manually determined)
# dendro.sign %>% plot(main='example', xlim=c(120,200), ylim=c(0,5))
# dendro.sign %>% plot(main='example', xlim=c(1300,1700), ylim=c(0,10))
# 
# for (i in 1:total){
#   if (i %in% signif[!is.na(signif)]){
#     if (i > (internal_nodes)){
#       node_colors[i] <- 'blue'
#     }
#     else{
#       children <- hc$merge[i,]
#       if (children[1]<0){
#         child1 <- internal_nodes + abs(children[1])
#       }
#       else{
#         child1 <- children[1]
#       }
#       if (children[2]<0){
#         child2 <- internal_nodes + abs(children[2])
#       }
#       else{
#         child2 <- children[2]
#       }
#       
#       ## Children aren't significant, combination is
#       if (R$pvals[child1] > 0.05 & R$pvals[child2] > 0.05){
#         node_colors[i] <- "purple"
#       }
#       ## Children are significant
#       if(R$pvals[child1] < 0.05 & R$pvals[child2] < 0.05){
#         
#         ## Parent contains more information by combining children
#         if (R$BIC[i] < R$BIC[child1] & R$BIC[i]<R$BIC[child2]){
#           node_colors[i] <- "purple"
#         }
#         
#         ## Children are more significant, node doesn't add anything (rare)
#         if(R$BIC[i] > R$BIC[child1] & R$BIC[i]>R$BIC[child2]){
#           node_colors[i] <- "orange"
#         }
#         
#         ##  Node is significant because of one child
#         if(R$BIC[i] < R$BIC[child1] | R$BIC[i]<R$BIC[child2]){
#           node_colors[i] <- "green"
#           print(i)
#         }
#       }
#       
#       ## One child is significant
#       if (R$pvals[child1] < 0.05 | R$pvals[child2] < 0.05){
#         sig_child <- if (R$pvals[child1] < R$pvals[child2]) child1 else child2
#         
#         if(R$BIC[i] < R$BIC[sig_child]) {
#           node_colors[i] <- "purple"
#         }
#         else{
#           node_colors[i] <- "orange"
#         }
#       }
#     }
#   }
# }
# 
# 
# platforms <- unique(rowData(SE)$platform)
# col_plats <- c("red", "yellow", "blue", "green", "purple", "orange", "darkgreen", "brown", "magenta", "pink", "grey")
# colors <- mapply(function(x) col_plats[which(platforms==x)], rowData(SE)$platform)
# colored_bars(colors, dendro.sign, add = T, rowLabels = c("Platform"))
# #### show module with lowest pvalue ----
# pdf(file = "purpl_ggm.pdf")
# for (node in purple){
#   plot.new()
#   members <- ff_2(hc$merge, node)
# 
#   # look at pearson correlation matrix of this cluster
#   # verify that these nodes are in one cluster (they must be adjacent in the hclust index list)
#   view_range <- match(members, hc$order)
#   par(fig = c(0,0.55,0,1),new = T, mar=c(2,3,2,.05))
#   dendro.sign %>% plot(main=node, xlim=c((min(view_range)-3),(max(view_range)+3)), ylim=c(0,1))
#   
#   
#   par(fig = c(0.5,1,0,1),new = T, mar=c(0,0,2,0))
#   adj_mat <-get_adj_mat(node, hc$merge, partial_correlation_matrix)[[1]]
#   g <- network(adj_mat)
#   plot.network(main =get_adj_mat(cor_vals)[[2]], g, usearrows = F, edge.lwd = 0.1, vertex.cex = 0.7)
# 
#   
# }
# dev.off()
# 
# get_adj_mat <- function(
#   node_number,
#   me,
#   correlation_mat
#   )
#   {
#   members <- ff_2(me, node_number)
#   cor_vals <- correlation_mat[members,members]
#   full_conn <- graph_from_adjacency_matrix((-1*abs(cor_vals)), mode='undirected', weighted = T)
#   min_span <- mst(full_conn, weights = E(full_conn)$weight)
#   cutoff <- abs(max(E(min_span)$weight))
#   adj_mat <- 1*(abs(cor_vals) >= cutoff)
#   adj_mat[lower.tri(adj_mat)] <- 0
#   list(adj_mat, cutoff)
# }
# 
# 
# 
