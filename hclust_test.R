library(magrittr)
library(survival)
library(dplyr)
library(ggplot2)

dat = t(assay(all_platforms))
# correlation to distance matrix
d = as.dist( 1-abs(cor(dat))) #exp(1-abs(cor(dat)) )-1 
# min threshold for clusters
th_min = 2
# max threshold for clusters
th_max = 50
  
re<- seq(0.15,1, 0.05) %>% # L-W alphas 
  lapply( function(linkage)try({
    # equivalently if alpha is a linkage i.e. "median", "single"
    # hc = hclust(d, method = linkage)
    
    hc = as.hclust(cluster::agnes(d, diss = T,method = "flexible", par.method = list(linkage)))
    # cut the tree for every level possible 
    khc = cutree(hc, k = 1:ncol(dat))[hc$order,]
    # some levels contain no clusters we might be interested in 
    khc= khc[ ,apply( khc ,2, function(x) min(table(x)) )< (th_max+1)]
    khc= khc[ ,apply( khc ,2, function(x) max(table(x)) )> (th_min)]
    
    
    # cluster labeller, put an unique label to each cluster, kind of hashing  
    f_clab <- function(v){ 
      sapply(sort(unique(v)),function(i){
        # v = sort(v)
        ins = which(v==i)
        #browser()
        if(length(ins)<th_min) return("")
        if(length(ins)>th_max) return("")
        paste(paste(min(ins), max(ins), sep = "_"), i, sep = ".")
      })
    }
    
    # get unique clusters 
    cls = apply(khc,2,f_clab) %>% unlist %>% unique
    # remove invalid hash ids 
    cls = cls[cls!=""]
    
    # cluster positions from hashed ids
    df = lapply( cls, strsplit, split="_|\\.") %>% 
      lapply(function(x) as.numeric(unlist(x))) %>% 
      do.call(what = rbind) 
    
    colnames(df) = c("i1", "i2", "id")
    df = data.frame(df) %>% arrange(i2-i1) 
    
    # remove if still small clusters there
    df = df[(df$i2-df$i1) > th_min,]
    df = df[,c("i1","i2")] %>% unique
    df = df %>% arrange(i1)
    
    # reduce the overlapping clusters, remove all the clusters and 
    # keep only largest subtree 
    
    f_red <- function(df){
      mm = as.matrix(df)
      mm = mm[order(mm[,2]-mm[,1], decreasing = T),]
      lind = lapply(seq(nrow(mm)), function(i) setdiff(which(((mm[,1]>= mm[i,1]) + ( mm[,2]<=mm[i,2]))==2),i))
      mm[-unique(unlist(lind)), ] 
    }
    df = data.frame(f_red(df))%>% arrange(i1)
    
    # pearson correlation between whole tree coph and distances 
    coph0 = cor( as.dist(cophenetic(hc)) %>% as.vector, as.vector(d))
    # kendall 
    cophx =  2*concordance( as.vector(as.dist(cophenetic(hc)))~as.vector(d))$concordance-1
    
    # measure for each subtree
    coph1 =
      lapply(seq(nrow(df)), function(k){
        ii = as.matrix(df)[k,]
        i = ii[1]:ii[2]
        # hc distances
        d_hc =cophenetic(hc) %>% as.matrix %>% {as.dist(.[i,i])} %>% as.vector
        # original distances
        d_ma = as.matrix(d)[i,i] %>%  as.dist %>% as.vector
        # cor(d_hc, d_ma)
        data.frame(d_hc, d_ma, k = k, l =length(d_hc))
      }) %>% do.call(what = rbind)
    
    # kendall for weighted mean of each subtree
    xkend = survival::concordance(d_hc~d_ma+strata(k), coph1)$concordance %>% {2*.-1}
    # pearson corrs of subtrees 
    xpear = coph1 %>% group_by(k) %>% summarise(xc = cor(d_hc, d_ma), n= unique(l))
    list(xkend = xkend, xpear =xpear, pear0 =coph0, kend0 = cophx  )
    
  }))

# drop if some setting does not work 
re = re[sapply(re, is.list)]

# get the summary from the results 
sm = sapply(re, function(x) c(unlist( x[c("xkend","pear0","kend0")] ), xpear = median(x$xpear$xc) ))

sm2 = lapply(re, function(x) x$xpear$xc)
names(sm2) = paste0(names(sm2), "_")
sm2  = unlist(sm2)
sm2 = data.frame(v=sm2, name = gsub("_.*","", names(sm2)))


mdf = reshape2::melt(sm)
# what kind of tree
mdf$Q <- c("subtree","alltree")[ 2-(mdf$Var1 %in% c("xkend", "xpear"))]
# what kind of corr coefficient
mdf$P <- c("kendall","pearson")[ 2-(mdf$Var1 %in% c("kend0", "xkend"))]

# library(ggplot2)
# ggplot(mdf, aes(x = Var2, y = value)) + 
#   geom_point() + 
#   geom_line(aes(x = as.numeric(Var2), lty = Q, group = Var1, color = P)) +
#   geom_boxplot(data = sm2, mapping = aes(x=name, y =v), width = 0.1, height = 0.1)

# sm2: individual subtrees
# mdf: summary results
ggs <- list(mdf = mdf, sm2 = sm2)


# combine all results into giant data frame 

colnames(ggs$mdf) <- c("setting", "alpha", "coph","tree", "method")
colnames(ggs$sm2) <- c("coph", "alpha")
    ggs$sm2$tree = "subtree"
    ggs$sm2$method = "pearson"
    ggs$sm2$setting = "single_tree"
    rbind(ggs$mdf, ggs$sm2[, colnames(ggs$mdf)])

# needed for ggplot
median_se<-function(x) {
  x <- stats::na.omit(x)
  bp = boxplot(x,plot = F, notch = T)[c("stats","conf")]
  nt = bp$conf
  st = bp$stats
  ggplot2:::new_data_frame(list(y = st[3], ymin = nt[1], ymax = nt[2]), n = 1)
}

# giant data frame
df <- ggs
df$alpha <- as.numeric(df$alpha)
df$method = factor(df$method, levels = c("pearson", "kendall"))
# which results are you interested in
df = df[df$method=="pearson",]
# plot
ggplot(df[df$setting!="single_tree",], aes(x = alpha, y = coph, color = tree)) +
  # geom_boxplot(data = df[df$setting=="single_tree",], aes(group = alpha),
  #              fill =NA, outlier.color = NA )+
  # stat_summary(data = df[df$setting=="single_tree",],aes(group = alpha),
  #              fun.y=mean, geom="point", size=2) + 
  geom_vline(data = df[df$setting!="single_tree",] %>% 
               group_by(cancer, tree) %>% 
               summarise(m =  alpha[which.max(coph)]),
             mapping = aes(xintercept = m, color = tree), lwd = 1.2, alpha = 0.5)+
  stat_summary(data = df[df$setting=="single_tree",],aes(group = alpha),
               fun.data = median_se, geom = "errorbar",width=0.01, 
               alpha = 0.25)+
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) + 
  facet_wrap(.~cancer,scales = "free_y",nrow = 3) + 
  theme_minimal() + 
  theme(strip.background = element_rect(fill = "wheat", color = "black"),
        panel.background = element_rect(colour = "black") ) + 
  ggsci::scale_color_aaas() +
  theme(legend.position = c(0.8,0.2))+
  scale_y_continuous(limits = c(0,NA))+
  labs(x= "Lance-Williams alpha", y = "pearson( cophenetic(hc), 1-abs(cor(data)) )",
       title = "Comparing cophenetic consistency of proposed subtree-version vs whole-tree")

# kendall plot
df = do.call(rbind, aa)
df$alpha <- as.numeric(df$alpha)
df$method = factor(df$method, levels = c("pearson", "kendall"))
df = df[df$method!="pearson",]
ggplot(df[df$setting!="single_tree",], aes(x = alpha, y = coph, color = tree)) +
  # geom_boxplot(data = df[df$setting=="single_tree",], aes(group = alpha),
  #              fill =NA, outlier.color = NA )+
  # stat_summary(data = df[df$setting=="single_tree",],aes(group = alpha),
  #              fun.y=mean, geom="point", size=2) + 
  geom_vline(data = df[df$setting!="single_tree",] %>% 
               group_by(cancer, tree) %>% 
               summarise(m =  alpha[which.max(coph)]),
             mapping = aes(xintercept = m, color = tree), lwd = 1.2, alpha = 0.5)+
  stat_summary(data = df[df$setting=="single_tree",],aes(group = alpha),
               fun.data = median_se, geom = "errorbar",width=0.01, 
               alpha = 0.25)+
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) + 
  facet_wrap(.~cancer,scales = "free_y",nrow = 3) + 
  theme_minimal() + 
  theme(strip.background = element_rect(fill = "wheat", color = "black"),
        panel.background = element_rect(colour = "black") ) + 
  ggsci::scale_color_aaas() +
  theme(legend.position = c(0.8,0.2))+
  scale_y_continuous(limits = c(0,NA))+
  labs(x= "Lance-Williams alpha", y = "kendall( cophenetic(hc), 1-abs(cor(data)) )",
       title = "Comparing cophenetic consistency of proposed subtree-version vs whole-tree")
