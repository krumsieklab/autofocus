## Testing AutoFocus##

suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(glmnet))

# Preprocess QMDiab data
source(codes.makepath("autofocus/PreprocessQMDiab.R"))

# Preprocess ADNI data
source(codes.makepath("autofocus/ADNI_data_loading.R"))
source(codes.makepath("autofocus/Backend_Modules.R"))
source(codes.makepath("autofocus/initialize.R"))
source(codes.makepath("autofocus/ModuleTesting.R"))

#set adjust_method 
#wy for westfall young - takes forever
#any other adjustment method that p.adjust takes
adjust_method = "holm"

# set number of cores
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


AD_results_pc <- find_sig_clusts(AD, 
                              "ADAS.Cog13", c("Age", "Sex", "Education"), 
                              score_method = "pc", 
                              adjust_method = adjust_method,
                              cores = cores)

