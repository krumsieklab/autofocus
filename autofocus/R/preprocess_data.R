#' Bind SE objects - handle NAs for non-overlapping samples
#'
#' Binds the SummarizedExperiments of individual platforms into a single SummarizedExperiment where
#' all samples are kept and NAs are (purposefully) produced where a sample doesn tn have an
#' analyte measurement.
#' Note that the sample names must be in the same format and under the same column name for all
#' platforms!
#'
#' @param platforms List of individual platform SEs.
#' @param sample_id Column name of the sample identifier.
#'
#' @return A SummarizedExperiment object with all platforms data
#'
#' @export
bind_se_with_na <- function(platforms, sample_id){

  # Set the column names of each platform to the universal sample identifier
  platforms <- lapply(platforms, function(i) {
    colnames(i)<-colData(i)[[sample_id]]
    i}
  )

  # Bind the assay data, order
  full_assay <- lapply(platforms, function(i) data.frame(assay(i), check.names=FALSE)) %>% bind_rows()

  columns_df <- lapply(platforms, function(i) data.frame(colData(i)))

  full_colData <- Reduce(function(x, y) merge(x, y, by=colnames(columns_df[[1]]), all=T), columns_df)
  rownames(full_colData) <- full_colData$SampleId
  full_colData <- full_colData[colnames(full_assay),]

  full_rowData <- matchRowClasses(platforms)

  SummarizedExperiment(assays = full_assay, colData = full_colData, rowData = full_rowData)
}

#' Combine rowData data frames from list of platforms
#'
#' Combine rowData data frames from each SE object in a list of platforms.
#'
#' @param platforms List of individual platform SEs.
#'
#' @return Combined rowData ata frame
#'
#' @noRd
matchRowClasses <-  function(platforms){

  rows_df <- lapply(platforms, function(i) data.frame(rowData(i)))

  row_names_types <- lapply(rows_df, function(x) {
    data.frame(Name = names(x), Type = sapply(x, typeof))
  }) %>%
    bind_rows()  %>%
    group_by(Name) %>%
    slice(1) %>%
    ungroup() %>%
    droplevels()

  for(x in 1:length(rows_df)){
    for (i in 1:ncol(rows_df[[x]])){
      class(rows_df[[x]][,i])<- as.character(row_names_types$Type[which(row_names_types$Name == colnames(rows_df[[x]])[i])])
    }
  }

  rows_df %>% bind_rows()

}
