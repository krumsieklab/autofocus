source(codes.makepath("autofocus/internal_funcs.R"))
source(codes.makepath("autofocus/hierarchical.R"))

set.seed(NULL)

make_R <- function(
  SE,
  mol_data_rows,
  phenotype,
  confounders,
  save_file = F,
  filename = "outfile.rmd",
  cores = 7
){
  data_mat <- t(assay(SE))
  sample_data <- data.frame(colData(SE))
  mol_data <- data.frame(rowData(SE))[mol_data_rows]
  R <- initialize(data_mat, sample_data, mol_data) %>% 
    find_sig_clusts(phenotype, confounders, cores)
  if(save_file){
    save(R, file = filename)
  }
  R
}


##################################################################


