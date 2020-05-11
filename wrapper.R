source(codes.makepath("autofocus/internal_funcs.R"))
source(codes.makepath("autofocus/hierarchical.R"))

set.seed(NULL)

make_R <- function(
  SE,
  mol_data_rows,
  phenotype,
  confounders,
  filename,
  cores = 7
){
  data_mat <- t(assay(SE))
  sample_data <- data.frame(colData(SE))
  mol_data <- data.frame(rowData(SE))[mol_data_rows]
  R <- initialize(data_mat, sample_data, mol_data) %>% 
    find_sig_clusts(phenotype, confounders, cores)
  save(R, file = filename)
  R
}



# If we want significance test to be performed

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


