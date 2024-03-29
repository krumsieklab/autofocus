#' Merge List of SE objects
#'
#' Merges a list of SE objects provided/ Combines assay, rowData, and colData data frames. Each
#' SE object must satisfy the following requirements:
#' \itemize{
#'    \item{the name of the colData column containing sample identifiers must be the same across all SE objects}
#'    \item{the name of the rowData column containing feature identifiers must be the same across all SE objects}
#'    \item{samples can overlap, but all features must be unique}
#'    \item{the column names in \code{phenos} and \code{confounders} must be present in the colData data frame of each SE object}
#' }
#'
#' @param D_list List of SE objects.
#' @param samp_id Name of colData column containing sample identifiers. Must be the same for each
#'    SE object.
#' @param feat_id Name of rowData column containing feature identifiers. Must be the same for each
#'    SE object.
#' @param phenos A vector containing the phenotypes of interest. These must be colData columns and
#'    they must be present in each SE object. This function will crash if these are not present.
#' @param confounders A vector containing confounders. These must be colData columns and they must
#'    be present in each SE object. This function will crash if these are not present.
#'
#' @return Merges list of SE objects into a single large SE object.
#'
#' @examples
#' \dontrun{D_all <-
#'   # merge list of SE objects
#'   merge_se(D_list = list(D1, D2),
#'   samp_id="Sample_Id",
#'   feat_id="names",
#'   phenos="DIAB",
#'   confounders=c("AGE", "SEX", "BMI")) %>%
#'   ...}
#'
#' @author KC
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment assay
#'
#' @export
merge_se <- function(D_list, samp_id, feat_id, phenos, confounders){

  # Check that samp_id and assay columns match; throw warning if not
  sample_list <- sapply(D_list, function(x){colData(x)[[samp_id]]}) %>% unlist()
  assay_cols <- sapply(D_list, function(x){colnames(assay(x))}) %>% unlist()
  if(!identical(sample_list, assay_cols)) warning("Assay column names do not match sample IDs. User must ensure correct ordering of assay columns.")

  # check features are unique
  feature_list <- sapply(D_list, function(x){rowData(x)[[feat_id]]}) %>% unlist()
  if(any(duplicated(feature_list))) stop("Features must be unique!")

  # Check that feat_id and assay rows match; throw warning if not
  assay_rows <- sapply(D_list, function(x){rownames(assay(x))}) %>% unlist()
  if(!identical(feature_list, assay_rows)) warning("Assay row names do not match feature IDs. User must ensure correct ordering of assay rows.")

  # check required colData columns are present
  req_col <- c(samp_id, phenos, confounders)
  check_cols <- lapply(D_list, function(x){req_col %in% colnames(colData(x))}) %>% unlist()
  if(!any(check_cols)) stop("Phenotype and / or confounder columns missing in one or more SEs.")

  # get combined assay
  df_list <- lapply(D_list, function(x){
    df <- assay(x) %>% as.data.frame()
    colnames(df) <- colData(x)[[samp_id]]
    rownames(df) <- rowData(x)[[feat_id]]
    df <- df %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(samp_id)
    df
  })
  full_df <- df_list %>% purrr::reduce(dplyr::full_join, by=samp_id) %>% tibble::column_to_rownames(samp_id) %>%
    t() %>% as.data.frame()

  # sample colData
  cd_list <- lapply(D_list, function(x){colData(x) %>% as.data.frame() %>% dplyr::select(dplyr::all_of(req_col))})
  full_cd <- data.table::rbindlist(cd_list) %>% dplyr::distinct()

  # feature rowData
  rd_list <- lapply(D_list, function(x){rowData(x) %>% as.data.frame()})
  full_rd <- data.table::rbindlist(rd_list, fill = TRUE)

  # create new SummarizedExperiment
  D <- SummarizedExperiment::SummarizedExperiment(assay=full_df,
                                                  colData=full_cd,
                                                  rowData=full_rd,
                                                  metadata = list())

  D

}
