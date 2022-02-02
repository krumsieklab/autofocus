#### Initialize ----
# zap()
# load libraries
library(MetaboTools)
library(tidyverse)
library(magrittr)
library(dplyr)
#### Define parameters ----
# max missingness for quotient normalization reference sample calculation
max_miss_norm <- 0.2
# max missingness for filtering
max_miss_filter <- 0.5
#### Load Data ----
file_data <- data.makepath("shareddata/ProstateCancer/EDRN/CORN-02-19ML+/CORN-02-19ML+ CDT_2019-12-23_MTfriendly.xlsx")
file_clin <- data.makepath("shareddata/ProstateCancer/EDRN/EDRN_ClinicalInfo_2020-01-10.xlsx")
# load data
D <- 
  #mt_reporting_heading(strtitle = "Load Data", lvl = 1) %>%
  # load data
  mt_files_data_xls(file=file_data, sheet="data", samples_in_row=T, ID_col="SAMPLE_ID") %>% 
  # load sample annotations
  mt_files_anno_xls(file=file_data, sheet="sampleinfo", anno_type="samples", anno_ID="SAMPLE_ID") %>% 
  # load metabolite annotations
  mt_files_anno_xls(file=file_data, sheet="metinfo",anno_type="metabolites", anno_ID="BIOCHEMICAL", data_ID="name") %>%
  # load clinical annotations
  mt_files_anno_xls(file=file_clin, sheet="clin", anno_type="samples", anno_ID="A1", data_ID="SAMPLE_ID") %>%
  {.}
#### Preprocessing ----
D <- D %>%
  mt_reporting_heading(strtitle = "Preprocessing", lvl = 1) %>%
  # filtering
  mt_plots_qc_missingness(met_max=max_miss_filter) %>%
  mt_pre_filtermiss(met_max=max_miss_filter) %>%
  mt_plots_qc_missingness(met_max=max_miss_filter) %>%
  # sample boxplots
  mt_plots_sampleboxplot(color=Diagnosis, plottitle = "Original", logged = T) %>%
  # normalization
  mt_pre_norm_quot(met_max = max_miss_norm) %>% 
  mt_plots_qc_dilutionplot(comp="Diagnosis") %>%
  # sample boxplots
  mt_plots_sampleboxplot(color=Diagnosis,plottitle = "After normalization", logged = T) %>%
  # log transformation
  mt_pre_trans_log() %>%
  # imputation
  mt_pre_impute_knn() %>%
  {.}
#### Add KEGG pathways ----
D %<>% 
  mt_anno_pathways_HMDB(in_col = "HMDB", out_col = "KEGG_pw", pwdb_name = "KEGG")
#### Extract elements ----
# metabolomics dataframe
data <- D %>% assay %>% t
# phenotype dataframe (Prostate Cancer Diagnosis 1=Yes, 0=No)
pheno <- D %>% colData %>% as.data.frame %>% dplyr::select(Diagnosis,Age)
# metabolite annotation dataframe
metanno <- D %>% rowData %>% as.data.frame %>% dplyr::select(name,SUB.PATHWAY,SUPER.PATHWAY, KEGG_pw)

adjust_method = "wy"

# set number of cores
cores = 40


EDRN<- initialize_R(data, 
                    sample.data = pheno,
                    mol.data = metanno)

EDRN_results <- find_sig_clusts(EDRN, 
                                "Diagnosis", c("Age"), 
                                score_method = "lm", 
                                adjust_method = adjust_method, 
                                cores = cores)
