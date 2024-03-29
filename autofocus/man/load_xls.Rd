% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_xlsx.R
\name{load_xls}
\alias{load_xls}
\title{Load data and annotaitons from Excel file and create SE object}
\usage{
load_xls(file, data_sheet, feature_sheet, sample_sheet)
}
\arguments{
\item{file}{Name of input Excel file.}

\item{data_sheet}{Name of sheet corresponding to assay data frame.}

\item{feature_sheet}{Name of sheet corresponding to rowData data frame.}

\item{sample_sheet}{Name of sheet corresponding to colData data frame.}
}
\description{
Load data and annotation matrices from Excel sheet. This function assumes the Excel file
contains three sheets corresponding to the assay, rowData, and colData data frames. The file
must satisfy the following requirements:
\itemize{
   \item{The data sheet must be in samples-as-columns format, with the first column of the
   data frame corresponding to the feature IDs and the first row corresponding to sample IDs.}
   \item{The feature sheet must contain at least two columns each with the feature IDs (the
    first column will be converted to rownames).}
    \item{The sample sheet must contain at least one column with sample IDs.}
}
}
\examples{
\dontrun{
   D <- load_xls(file = "process_BM.xlsx",
                 data_sheet = "data",
                 feature_sheet = "feature_anno",
                 sample_sheet = "sample_anno")
}

}
\author{
KC
}
