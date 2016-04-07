#!/usr/bin/env Rscript
library(rmarkdown)
args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("Please specify the path to .Rmd file", call.=FALSE)
}
rmd_fpath = args[1]
sink(stdout(), type = "message")
sink("/dev/null")
render(rmd_fpath)
sink()