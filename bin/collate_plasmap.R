#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(impoRt)
library(tidyr)

data <- get_data(".",
                 "report.txt",
                 recursive = FALSE) %>%
  select(-ref) %>%
  rowwise() %>%
  mutate(id = sub(paste0(plasmid, "_"), "", id),
         id = sub("_mapped_sorted", "", id),
         percent_mapped = round(percent_mapped, 3)) %>%
  ungroup() %>%
  spread(plasmid, percent_mapped)

write.table(data,
            "plasmid_mapping_report.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

