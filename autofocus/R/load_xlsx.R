#' Load data and annotaitons from Excel file and create SE object
#'
#' Load data and annotation matrices from Excel sheet. This function assumes the Excel file
#' contains three sheets corresponding to the assay, rowData, and colData data frames. The file
#' must satisfy the following requirements:
#' \itemize{
#'    \item{The data sheet must be in samples-as-columns format, with the first column of the
#'    data frame corresponding to the feature IDs and the first row corresponding to sample IDs.}
#'    \item{The feature sheet must contain at least two columns each with the feature IDs (the
#'     first column will be converted to rownames).}
#'     \item{The sample sheet must contain at least one column with sample IDs.}
#' }
#'
#' @param file Name of input Excel file.
#' @param data_sheet Name of sheet corresponding to assay data frame.
#' @param feature_sheet Name of sheet corresponding to rowData data frame.
#' @param sample_sheet Name of sheet corresponding to colData data frame.
#'
#' @examples
#' \dontrun{
#'    D <- load_xls(file = "process_BM.xlsx",
#'                  data_sheet = "data",
#'                  feature_sheet = "feature_anno",
#'                  sample_sheet = "sample_anno")
#' }
#'
#' @author KC
#'
#' @export
load_xls <- function(file,
                      data_sheet,
                      feature_sheet,
                      sample_sheet){

  ### validate arguments ------
  if(missing(file)) stop("File must be provided.")
  if(missing(data_sheet)) stop("A name must be provided for argument \'data_sheet\'.")
  if(missing(feature_sheet)) stop("A name must be provided for argument \'feature_sheet\'.")
  if(missing(sample_sheet)) stop("A name must be provided for argument \'sample_sheet\'.")

  ### get assay, rowData, and colData matrices ------
  df <- as.data.frame(readxl::read_excel(path=file,sheet=data_sheet, col_names=T))
  df %<>% tibble::column_to_rownames(colnames(df)[1])
  assay <- as.matrix(df)
  rd <- as.data.frame(readxl::read_excel(path=file,sheet=feature_sheet,col_names=T))
  rd %<>% tibble::column_to_rownames(colnames(rd)[1])
  cd <- as.data.frame(readxl::read_excel(path=file,sheet=sample_sheet,col_names=T))

  ### construct SummarizedExperiment ------
  D <- SummarizedExperiment(assay = assay,
                            rowData = rd,
                            colData = cd,
                            metadata = list(sessionInfo=utils::sessionInfo()))

  D

}
