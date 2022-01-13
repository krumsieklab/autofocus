#rm(list=ls())
# alzheimer data
library(dplyr)
# reading data sets

phe =  read.csv2(data.makepath("shareddata/ADNI/ADNI_datasets_from_Matthias/adni1go2.phenotypes.covariates.csv"), na.string="NA",header=T,sep=",",dec = ".")
bba =  read.csv2(data.makepath("shareddata/ADNI/ADNI_datasets_from_Matthias/adni1go2.ba.medadjusted.csv"), na.string="NA",header=T,sep=",",dec = ".")
bp1 =  read.csv2(data.makepath("shareddata/ADNI/ADNI_datasets_from_Matthias/adni1go2.p180.medadjusted.csv"), na.string="NA",header=T,sep=",",dec = ".")
ntl =  read.csv2(data.makepath("shareddata/ADNI/ADNI_datasets_from_Matthias/adni1.lipids.medadjusted.csv"), na.string="NA",header=T,sep=",",dec = ".")
mkl = read.csv2(data.makepath("shareddata/ADNI/ADNI_datasets_from_Matthias/adni1.meikle.medadjusted.csv"), na.string="NA",header=T,sep=",",dec = ".")


bind_SE_no_NA <- function(
  platforms  # SummarizedExperiments to combine
){
  for (platform in platforms){
    if (exists("dfTot")){
      missing1 <- setdiff(colnames(rowData(dfTot)), colnames(rowData(platform)))
      missing2 <- setdiff(colnames(rowData(platform)), colnames(rowData(dfTot)))
      for (colname in missing1){
        rowData(platform)[[colname]] <- rep(NA, nrow(platform))
      }
      for (colname in missing2){
        rowData(dfTot)[[colname]]<- rep(NA, nrow(dfTot))
      }
      
      overlap <- intersect(colnames(dfTot), colnames(platform))
      dfTot <-  rbind(dfTot[,overlap], platform[,overlap])
    }
    else{
      dfTot <-platform
    }
  }
  dfTot
}

bba.overlap <- intersect(phe$RID, bba$RID)
SE.bba <- SummarizedExperiment(assays = t(bba[,3:ncol(bba)]), 
                               colData = phe[which(phe$RID %in% bba.overlap),])
colnames(SE.bba)<-bba[,1]
rowData(SE.bba)$name <- colnames(bba)[3:ncol(bba)]
rowData(SE.bba)$platform <- rep('bba', times = (ncol(bba)-2))

bp1.overlap <- intersect(phe$RID, bp1$RID)
SE.bp1 <- SummarizedExperiment(assays = t(bp1[,3:ncol(bp1)]), 
                               colData = phe[which(phe$RID %in% bp1.overlap),])
colnames(SE.bp1)<-bp1[,1]
rowData(SE.bp1)$name <- colnames(bp1)[3:ncol(bp1)]
rowData(SE.bp1)$platform <- rep('bp1', times = (ncol(bp1)-2))

ntl.overlap <- intersect(phe$RID, ntl$RID)
SE.ntl <- SummarizedExperiment(assays = t(ntl[,2:ncol(ntl)]), 
                               colData = phe[which(phe$RID %in% ntl.overlap),])
colnames(SE.ntl)<-ntl[,1]
rowData(SE.ntl)$name <- colnames(ntl)[2:ncol(ntl)]
rowData(SE.ntl)$platform <- rep('ntl', times = (ncol(ntl)-1))

mkl.overlap <- intersect(phe$RID, mkl$RID)
SE.mkl <- SummarizedExperiment(assays = t(mkl[,2:ncol(mkl)]), 
                               colData = phe[which(phe$RID %in% mkl.overlap),])
colnames(SE.mkl)<-mkl[,1]
rowData(SE.mkl)$name <- colnames(mkl)[2:ncol(mkl)]
rowData(SE.mkl)$platform <- rep('mkl', times = (ncol(mkl)-1))

ADNI_platforms <- bind_SE_no_NA(list(SE.bba, SE.bp1, SE.ntl, SE.mkl))%>% mt_pre_confounding_correction(formula = ~Sex+Age+bmi_in_kg_p_m2+Education)
