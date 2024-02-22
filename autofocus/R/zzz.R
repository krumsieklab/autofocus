.onLoad <- function(libname, pkgname) {
  invisible(suppressPackageStartupMessages(
    sapply(c("cluster", "dendextend", "dplyr", "foreach", "ggplot2", "glmnet", "magrittr", "mgm", "parallel", "RColorBrewer", "reshape2", "SummarizedExperiment"),
           requireNamespace, quietly = TRUE)
  ))
}
