#### Initialize ----
# zap()
# load libraries
library(MetaboTools)
library(tidyverse)
library(magrittr)
library(dplyr)
source(codes.makepath("autofocus/Backend_Modules.R"))
source(codes.makepath("autofocus/initialize.R"))
source(codes.makepath("autofocus/ModuleTesting.R"))


####ADNI####

get_adni_R<- function(){
  
  # Load ADNI
  load(data.makepath("ADNI/ADNI_data_for_Annalise/ADNI_data_for_Annalise.RData"))
  source(codes.makepath("autofocus/ADNI_data_loading.R"))

  adni.data<- t(assay(ADNI_platforms))
  
  chemical_classes<- lapply(colnames(adni.data), function(i){names(fset.list)[grep(i,fset.list)[1]]}) %>% unlist()
  adni.metanno<-data.frame(rowData(ADNI_platforms), chemical_class=chemical_classes)
  
  
  ADNI <- initialize_R(adni.data,
                       sample.data = colData(ADNI_platforms),
                       mol.data = adni.metanno)
  ADNI
}

test_ADNI <- function(adni_R,score_method="pc",adjust_method="wy",cores=40){
  find_sig_clusts(adni_R,
                  "ADAS.Cog13", 
                  c("Age","Sex","bmi_in_kg_p_m2","Education"),
                  score_method,
                  adjust_method,
                  cores)
}



#### QMDIAB ####

get_qd_R<-function(){
  
  # Load QMDiab
  load(data.makepath("QMDiab/qmdiab_2019_03_13.rda"))
  source(codes.makepath("autofocus/PreprocessQMDiab.R"))
  
  names(all_platforms)<-rowData(all_platforms)$name
  QD <- initialize_R(t(assay(all_platforms)), 
                   sample.data = colData(all_platforms),
                   mol.data = rowData(all_platforms))
  QD
}

test_QMDiab <- function(qd_R, score_method="pc",adjust_method="wy",cores=40){
  find_sig_clusts(qd_R, 
                 "DIAB", 
                 c("AGE", "SEX", "BMI"), 
                 score_method, 
                 adjust_method, 
                 cores)
}



#### ROSMAP ####

get_rosmap_R<-function(){
  
  # Load ROSMAP
  load(data.makepath("ROSMAP/processed_data/autofocus/2021-03-04_ROSMAP_514Brains_AC.rds"))
  

  names(AC$D)<-rowData(AC$D)$BIOCHEMICAL

  ROSMAP<- initialize_R(t(assay(AC$D)), 
                      sample.data = colData(AC$D),
                      mol.data = rowData(AC$D))
  ROSMAP
}

test_ROSMAP<- function(rosmap_R,score_method="pc",adjust_method="wy",cores=40){
  find_sig_clusts(rosmap_R, 
                   "cogng_random_slope", 
                   score_method, 
                   adjust_method,
                   cores)
  }


