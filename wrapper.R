#source("PreprocessQMDiab.R")

source("internal_funcs.R")
source("hierarchical.R")
limit = 400000000
cores = 4
set.seed(NULL)
# input
### ADNI
SE <- ADNI_platforms
data_mat <- t(assay(SE))
sample_data <- data.frame(colData(SE))
mol_data <- data.frame(rowData(SE))[c(2)]
phenotype1 = "WholeBrain"
phenotype2 = "ADAS.Cog13"
confounders = c("Age","Sex","Education")
###

R <- initialize(data_mat, sample_data, mol_data) %>% 
  find_sig_clusts(R, phenotype1, confounders, cores)
save(R, file = "ADNI_wholebrain.rmd")

R <- initialize(data_mat, sample_data, mol_data) %>% 
  find_sig_clusts(R, phenotype2, confounders, cores)
save(R, file = "ADNI_adascog13.rmd")

###QMDIAB
SE <- bind_SE_no_NA(platform_list[5])
data_mat <- t(assay(SE))
sample_data <- data.frame(colData(SE))
mol_data <- data.frame(rowData(SE))[c(2,3,7,8)]
phenotype = "DIAB"
confounders = c("SEX","AGE","BMI")
###

R <- initialize(data_mat, sample_data, mol_data) %>% 
  find_sig_clusts(R, phenotype, confounders, cores)
save(R, file = "qmdiab_diab.rmd")

# workflow (outside code, "user view")
R <- initialize(data_mat, sample_data, mol_data) # creates a new R structure with $HCL, $clust_list, $C

# If we want significance test to be performed
R <- find_sig_clusts(R, phenotype, confounders, cores)
# after this, only R is ever needed, every function should add intermediate results to R

# s <- comp_size_of_cache(R)
# if (s <= limit) {
#   # size ok, we can cache
#   R <- R %>% 
#     cache() %>% 
#     TVM(FUN=cluster_net_cached) # make video
# } else {
#   # can't make the cache
#   R <- TVM(R, FUN=cluster_net)
# }


##################################################################


