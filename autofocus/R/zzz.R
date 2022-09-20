.onLoad <- function(libname, pkgname) {
  invisible(suppressPackageStartupMessages(
    sapply(c("cluster", "dendextend", "foreach", "glmnet", "magrittr", "mgm", "parallel", "SummarizedExperiment"),
           requireNamespace, quietly = TRUE)
  ))
}
