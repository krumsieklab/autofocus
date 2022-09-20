#### Initialize ----
# zap()
# load libraries
library(maplet)
library(tidyverse)
library(magrittr)
library(dplyr)
source(codes.makepath("autofocus/Backend_Modules.R"))
source(codes.makepath("autofocus/initialize.R"))
#source(codes.makepath("autofocus/ModuleTesting.R"))


####ADNI####

get_adni_R<- function(){
  
  # Load ADNI
  load(data.makepath("ADNI/ADNI_data_for_Annalise/ADNI_data_for_Annalise.RData"))
  
  # @KC: What is ADNI_platforms here?
  ADNI <- initialize_R(t(assay(ADNI_platforms)),
                       sample.data = colData(ADNI_platforms),
                       mol.data = rowData(ADNI_platforms),
                       confounders = c("Age","Sex","bmi_in_kg_p_m2","Education"),
                       phenotype="ADAS.Cog13")
  ADNI
}

dfs <- get_adni_R()
data.mat <- dfs$bp1[3:ncol(dfs$bp1)]
pheno.mat <- dfs$phe[dfs$phe$RID%in%dfs$bp1$RID,]
mol.mat <-  lapply(colnames(data.mat), function(i){names(fset.list)[grep(i,fset.list)[1]]}) %>% unlist()
BP1 <- initialize_R(data.mat[!is.na(pheno.mat$ADAS.Cog13),], mol.data=mol.mat, sample.data=pheno.mat[!is.na(pheno.mat$ADAS.Cog13),], confounders = c("Age", "Sex","Education","ApoE4","bmi_in_kg_p_m2"), "ADAS.Cog13")


#### QMDIAB ####

get_qd_R<-function(){
  
  # Load QMDiab
  load("/Users/aschweickart/Library/CloudStorage/Box-Box/shareddata/QMDiab/qmdiab_2019_03_13.rda")
  source(codes.makepath("autofocus/PreprocessQMDiab.R"))
  
  names(all_platforms)<-rowData(all_platforms)$name
  QD <- initialize_R(t(assay(all_platforms)), 
                   sample.data = colData(all_platforms),
                   mol.data = rowData(all_platforms),
                   c("AGE","SEX","BMI"),
                   "DIAB")
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
  load(data.makepath("shareddata/ROSMAP/processed_data/autofocus/rosmap_brain_3omics.rds"))
  D <- bind_SE_with_NA(list(Dm,Dp,Dt))
  medication_columns <- 'medication_columns.xlsx'
  Dmc <- D %>% # converting all medication columns to factors
    mt_anno_class(anno_type = 'samples', file=medication_columns, sheet=1) %>%
    # medication correction excluding AD related drugs
    mt_pre_confounding_correction_stepaic(cols_to_correct = names(colData(D))[grep("_rx", names(colData(D)))], cols_to_exclude = c("ad_rx", "neurologic_rx"))
  # confounder correction ----
  Dmc %<>% mt_anno_mutate(anno_type = "samples", col_name = 'msex',
                          term = as.factor(msex)) %>%
    mt_anno_mutate(anno_type = "samples", col_name = 'apoe_genotype',
                   term = as.factor(apoe_genotype)) %>%
    mt_anno_mutate(anno_type = "samples", col_name = 'pmi',
                   term = as.numeric(as.matrix(pmi))) %>%
    mt_anno_mutate(anno_type = "samples", col_name = 'age_death',
                   term = as.numeric(as.matrix(age_death))) %>%
    mt_anno_mutate(anno_type = "samples", col_name = 'bmi',
                   term = as.numeric(as.matrix(bmi))) %>%
    mt_anno_mutate(anno_type = "samples", col_name = 'educ',
                   term = as.numeric(as.matrix(educ)))

  ROSMAP<- initialize_R(t(assay(Dmc)), 
                      sample.data = colData(Dmc),
                      mol.data = rowData(Dmc),
                      confounders = c("msex","apoe_genotype","pmi","age_death","bmi","educ"),
                      "cogng_random_slope")
  ROSMAP
}

test_ROSMAP<- function(rosmap_R,score_method="pc",adjust_method="wy",cores=40){
  find_sig_clusts(rosmap_R, 
                   "cogng_random_slope", 
                   score_method, 
                   adjust_method,
                   cores)
  }


