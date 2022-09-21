library(dplyr)
library(maplet)
source(codes.makepath("autofocus/initialize.R"))
load(data.makepath("ROSMAP/processed_data/autofocus/rosmap_brain_3omics.rds"))
medication_columns <- 'medication_columns.xlsx'
names(Dm)<-rowData(Dm)$...1
names(Dp)<-rowData(Dp)$uname

Dmc <- Dm %>% # converting all medication columns to factors
  maplet::mt_anno_class(anno_type = 'samples', file=medication_columns, sheet=1) %>%
  # medication correction excluding AD related drugs
  maplet::mt_pre_confounding_correction_stepaic(cols_to_correct = names(colData(Dm))[grep("_rx", names(colData(Dm)))], cols_to_exclude = c("ad_rx", "neurologic_rx"))

Dpc <- Dp %>% # converting all medication columns to factors
  maplet::mt_anno_class(anno_type = 'samples', file=medication_columns, sheet=1) %>%
  # medication correction excluding AD related drugs
  maplet::mt_pre_confounding_correction_stepaic(cols_to_correct = names(colData(Dp))[grep("_rx", names(colData(Dp)))], cols_to_exclude = c("ad_rx", "neurologic_rx"))


platforms = list(Dmc,Dpc)
sample_id = "SAMPLE.NAME"
# colnames(Dm)<-colData(Dm)$SAMPLE.NAME
# colnames(Dp)<-colData(Dp)$SAMPLE.NAME
# full_assay<- lapply(platforms, function(i) data.frame(assay(i))) %>% bind_rows()
# full_assay<-full_assay[,order(as.numeric(colnames(full_assay)))]
# 
# columns_df <- lapply(platforms, function(i) data.frame(colData(i)))
# full_colData <- Reduce(function(x, y) merge(x, y, by=colnames(columns_df[[1]]), all=T), columns_df)
# full_colData <- full_colData[order(as.numeric(full_colData[[sample_id]])),]
# 
# 
# full_rowData <- lapply(platforms, function(i) data.frame(rowData(i))) %>% bind_rows()
# 
# D<- SummarizedExperiment(assays = full_assay, colData = full_colData, rowData = full_rowData)

D <- bind_SE_with_NA(platforms, sample_id)

ROSMAP<- initialize_R(t(assay(D)), 
                      sample.data = colData(D),
                      mol.data = rowData(D),
                      confounders = c("msex","apoe_genotype","pmi","age_death","bmi","educ"),
                      "cogng_random_slope")
