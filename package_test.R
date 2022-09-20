library(autofocus)
library(magrittr)
library(SummarizedExperiment)

load(data.makepath("ADNI/ADNI_data_for_Annalise/ADNI_data_for_Annalise.RData"))

data.mat <- dfs$bp1[3:ncol(dfs$bp1)]
pheno.mat <- dfs$phe[dfs$phe$RID%in%dfs$bp1$RID,]
mol.mat <-  lapply(colnames(data.mat), function(i){names(fset.list)[grep(i,fset.list)[1]]}) %>% unlist()
BP1 <- initialize_R(data.mat[!is.na(pheno.mat$ADAS.Cog13),], mol.data=mol.mat, sample.data=pheno.mat[!is.na(pheno.mat$ADAS.Cog13),], confounders = c("Age", "Sex","Education","ApoE4","bmi_in_kg_p_m2"), "ADAS.Cog13")


# QMDiab - urine or saliva dataset
# QMDiab or ROSMAP
glycan <- readRDS("/Users/kelsey/Downloads/glycan_data.rds")
confounders <- c("AGE", "SEX", "BMI")
data.mat <- assay(glycan)
pheno.mat <- data.frame(colData(glycan))
mol.mat <- data.frame(rowData(glycan))
R <- initialize_R(data.mat, pheno.mat, mol.mat, confounders = confounders)
