#!/usr/bin/env Rscript
library(dplyr)
library(funtools)
library(impoRt)

calc_cov <- function(df) {
  df %>%
    rename("reference" = X1,
           "pos" = X2,
           "reads" = X3) %>%
    mutate(reads_no_unmapped = ifelse(reads == 0, NA, reads),
           test = ifelse(reads == 0, 0, 1)) %>%
    summarise(reference = unique(reference),
              cov_perc = round(sum(test)/max(pos)*100, 2),
              mean_reads_per_base = round(mean(reads_no_unmapped, na.rm = TRUE), 2),
              median_reads_per_base = median(reads_no_unmapped, na.rm = TRUE),
              sd_reads_per_base = round(sd(reads_no_unmapped, na.rm = TRUE), 2))
}

cov_data <- get_data(".",
                     pattern = "genomecov.txt",
                     col_names = FALSE,
                     df = FALSE) %>%
  lapply(., calc_cov) %>%
  bind_rows(.id = "sample") %>%
  mutate(sample = sub("_genomecov.txt", "", sample))

write.table(cov_data,
            "plasmid_mapping_report.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)