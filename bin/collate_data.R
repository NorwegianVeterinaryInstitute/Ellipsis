#!/usr/bin/env Rscript
# ABSTRACT: Collate data from the plasmid pipeline into one table
library(dplyr)
library(impoRt)
library(funtools)

mobtyper_reports <- get_data(".",
                             pattern = "mobtyper",
			     convert = TRUE) %>%
  mutate(ref = sub("_mobtyper", "", ref),
         ref = sub(".fasta_report.txt", "", ref)) %>%
  select(-file_id)

plasmidfinder_reports <- get_data(".",
                                  pattern = "plasfinder",
				  convert = TRUE) %>%
  mutate(ref = sub("_plasfinder_results_tab.tsv", "", ref)) %>%
  group_by(ref) %>%
  summarise_all(list(func_paste)) %>%
  select(ref, Plasmid, Identity) %>%
  rename("pf_replicon" = Plasmid,
         "pf_replicon_identity" = Identity)

resfinder_reports <- get_data(".",
                              pattern = "resfinder",
                              convert = TRUE) %>%
  mutate(ref = sub("_resfinder_results_tab.tsv", "", ref)) %>%
  group_by(ref) %>%
  summarise_all(list(func_paste)) %>%
  select(ref, `Resistance gene`, Identity) %>%
  rename("rf_gene" = `Resistance gene`,
         "rf_gene_identity" = Identity)

virfinder_reports <- get_data(".",
                              pattern = "virfinder",
                              convert = TRUE) %>%
  mutate(ref = sub("_virfinder_results_tab.tsv", "", ref)) %>%
  group_by(ref) %>%
  summarise_all(list(func_paste)) %>%
  select(ref, `Virulence factor`, Identity) %>%
  rename("vf_gene" = `Virulence factor`,
         "vf_gene_identity" = Identity)

total_data <- mobtyper_reports %>%
  left_join(plasmidfinder_reports, by = "ref") %>%
  left_join(resfinder_reports, by = "ref") %>%
  left_join(virfinder_reports, by = "ref")

write.table(total_data,
            "total_report.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

