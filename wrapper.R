source(codes.makepath("autofocus/internal_funcs.R"))
source(codes.makepath("autofocus/ModuleTesting.R"))
source(codes.makepath("autofocus/app.R"))

args = commandArgs(trailingOnly = T)

if(length(args)<3){
  print("Not enough arguments supplied. Please include file, 
        --phenotypes tag followed by phenotype for testing,
         and --confounders tag followed by names of confounders")
}else{
  SE = load(args[1])
  phen_loc = which(args == "--phenotypes")
  if ("--confounders" %in% args){
    confounder_loc = which(args == "--confounders")
    phenotypes = args[(phen_loc+1):confounder_loc]
    confounders = args[(confounder_loc+1):length(args)]
  }
  else{
    phenotypes = args[(phen_loc+1):length(args)]
    confounders = c()
    }
}

if (length(phenotypes)==1){
  R <- make_R(SE, phenotypes, confounders)
}else{
  R_list <- mapply(phenotypes, function(p) make_R(SE, p, confounders))
}

run_App()


### Wrapper function taking in raw data, outputting significantly labelled dendrogram  ----------------

make_R <- function(
  SE,
  phenotype,
  confounders,
  save_file = F,
  filename = "outfile.rmd",
  cores = 7
){
  data_mat <- t(assay(SE))
  sample_data <- data.frame(colData(SE))
  mol_data <- data.frame(rowData(SE))
  R <- initialize(data_mat, sample_data, mol_data) %>% 
    find_sig_clusts(phenotype, confounders, cores)
  if(save_file){
    save(R, file = filename)
  }
  R
}

##################################################################


