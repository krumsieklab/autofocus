library(dplyr)
load(data.makepath("ROSMAP/processed_data/autofocus/rosmap_brain_3omics.rds"))
medication_columns <- 'medication_columns.xlsx'
colnames(Dm)<-colData(Dm)$SAMPLE.NAME
colnames(Dp)<-colData(Dp)$SAMPLE.NAME
platforms = list(Dm,Dp)
sample_id = "SAMPLE.NAME"
mol_id = "name"
full_assay <- lapply(platforms, function(i) data.frame(assay(i))) %>% bind_rows()
full_assay<-full_assay[,order(as.numeric(colnames(full_assay)))]

columns_df <- lapply(platforms, function(i) data.frame(colData(i)))
full_colData <- Reduce(function(x, y) merge(x, y, by=colnames(columns_df[[1]]), all=T), columns_df)
full_colData <- full_colData[order(as.numeric(full_colData[[sample_id]])),]

rows_df <- lapply(platforms, function(i) data.frame(rowData(i)))
full_rowData <- Reduce(function(x, y) merge(x, y, by=intersect(colnames(x),colnames(y)), all=T), rows_df)
full_rowData <- full_rowData[match(rownames(full_assay), full_rowData[[mol_id]]),]

D<- SummarizedExperiment(assays = full_assay, colData = full_colData, rowData = full_rowData)

Dmc <- D %>% # converting all medication columns to factors
  maplet::mt_anno_class(anno_type = 'samples', file=medication_columns, sheet=1) %>%
  # medication correction excluding AD related drugs
  maplet::mt_pre_confounding_correction_stepaic(cols_to_correct = names(colData(D))[grep("_rx", names(colData(D)))], cols_to_exclude = c("ad_rx", "neurologic_rx"))
