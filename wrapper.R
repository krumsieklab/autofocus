source(codes.makepath("autofocus/internal_funcs.R"))
source(codes.makepath("autofocus/ModuleTesting.R"))

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
    phenotypes = args[(phen_loc+1):(confounder_loc-1)]
    confounders = args[(confounder_loc+1):length(args)]
  }
  else{
    phenotypes = args[(phen_loc+1):length(args)]
    confounders = c()
    }
}

if (length(phenotypes)==1){
  R <- make_R(get(SE), phenotypes, confounders)
}else{
  R_list <- mapply(function(p) make_R(get(SE), p, confounders), phenotypes)
}

source(codes.makepath("autofocus/app.R"))




