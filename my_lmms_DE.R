my_lmmsDE<-function (data, time, sampleID, group, knots, numCores) 
{

 raw_ps <- mclapply(1:ncol(data), function(i){
   p_i <- try(get_model_ps(i, data, sampleID, time, group, knots))
   if (class(p_i) != "try-error") return(p_i) else return(NA)
    
  },mc.cores=numCores) %>% unlist()
  
  data.frame(Molecule=colnames(data), P_unadj=raw_ps, P_adj = p.adjust(raw_ps, method="fdr", n=length(raw_ps[!is.na(raw_ps)])))

}

get_model_ps<-function(i, data, sampleID, time, group, knots){
  dataM <- data.frame(Rep = sampleID)
  dataM$Group <- as.factor(group)
  dataM$time = as.numeric(time)
  dataM$Expr = as.numeric(data[,i])
  dataM$all = rep(1, nrow(dataM))

  PZ <- outer(dataM$time, knots, "-")
  PZ <- PZ * (PZ > 0)
  dataM$Zt <- PZ

  fit0L2 <- lme(Expr ~ as.factor(Group) + time, 
                data = dataM, 
                random = list(all = pdIdent(~Zt -1), 
                              Rep = pdDiag(~time)), 
                na.action = na.exclude, 
                control = lmeControl(opt = "optim"), 
                method = "ML")
  
  fit5L2 <- lme(Expr ~ as.factor(Group) + time + as.factor(Group) * time, 
                data = dataM, 
                random = list(Group = pdIdent(~Zt - 1), 
                              Rep = pdDiag(~time)), 
                na.action = na.exclude, 
                control = lmeControl(opt = "optim"), 
                method = "ML")

  
  anova.lme(fit0L2, fit5L2)$`p-value`[2][1]
}
