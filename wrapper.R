source(codes.makepath("autofocus/internal_funcs.R"))
source(codes.makepath("autofocus/ModuleTesting.R"))

# args = commandArgs(trailingOnly = T)
# if(length(args)<3){
#   print("Not enough arguments supplied. Please include file, 
#         --phenotypes tag followed by phenotype for testing,
#          and --confounders tag followed by names of confounders")
# }else{
#   SE = load(args[1])
#   phen_loc = which(args == "--phenotypes")
#   if ("--confounders" %in% args){
#     confounder_loc = which(args == "--confounders")
#     phenotypes = args[(phen_loc+1):(confounder_loc-1)]
#     confounders = args[(confounder_loc+1):length(args)]
#   }
#   else{
#     phenotypes = args[(phen_loc+1):length(args)]
#     confounders = c()
#     }
# }

# SE = ADNI_platforms
# phenotypes = c("WholeBrain", "ADAS.Cog13")
# confounders = c("Sex","Age","Education")
# 
# node_colors=c('red','darkorange','purple','green','blue')
# node_colors_light=c('pink','orange','mediumpurple','lightgreen','lightblue')
# 
# if (length(phenotypes)==1){
#   R <- make_R(SE, phenotypes, confounders, 'red', 'pink', save_file = T, filename = paste(phenotype,".rmd"))
# }else{
#   for (phenotype in phenotypes){
#     
#   }
#   R_list <- mapply(function(i) make_R(SE, phenotypes[1], confounders, node_colors[i], node_colors_light[i]), 
#                    1:length(phenotypes))
#   save(R_list, file = "results_list.rmd")
# }





