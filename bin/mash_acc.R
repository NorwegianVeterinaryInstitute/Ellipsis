#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
report_loc <- args[1]

library(readr)

report <- read_delim(report_loc, delim = "\t")

mash_acc <- report$mash_nearest_neighbor

print(mash_acc)

