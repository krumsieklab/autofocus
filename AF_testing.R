## Testing AutoFocus##

suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(epicalc))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(glmnet))

# source("PreprocessQMDiab.R")
# source("ADNI_data_loading.R")
# source("Backend_Modules.R")
# source("internal_funcs.R")
# source("ModuleTesting.R")



score_method = "lm"
#score_method = "pc"

adjust_method = "holm"
#adjust_method = "wy"

cores = 6


QD <- initialize(t(assay(all_platforms)), 
                   sample.data = colData(all_platforms),
                   mol.data = rowData(all_platforms))

QD_results_lm <- find_sig_clusts(QD, 
                              "DIAB", c("AGE", "SEX", "BMI"), 
                              score_method = "lm", 
                              adjust_method = adjust_method, 
                              cores = cores)

QD_results_pc <- find_sig_clusts(QD, 
                                 "DIAB", c("AGE", "SEX", "BMI"), 
                                 score_method = "pc", 
                                 adjust_method = adjust_method, 
                                 cores = cores)

AD <- initialize(t(assay(ADNI_platforms)), 
                 sample.data = colData(ADNI_platforms),
                 mol.data = rowData(ADNI_platforms))

AD_results_lm <- find_sig_clusts(AD, 
                              "ADAS.Cog13", c("Age", "Sex", "Education"), 
                              score_method = "lm", 
                              adjust_method = adjust_method,
                              cores = cores)


AD_results <- find_sig_clusts(AD, 
                              "ADAS.Cog13", c("Age", "Sex", "Education"), 
                              score_method = "pc", 
                              adjust_method = adjust_method,
                              cores = cores)

