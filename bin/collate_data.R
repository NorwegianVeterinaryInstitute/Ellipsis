#!/usr/bin/env Rscript
# ABSTRACT: Collate data from the plasmid pipeline into one table
library(dplyr)
library(impoRt)
library(funtools)
library(tidyr)

mobtyper_reports <- get_data(".",
                             pattern = "mobtyper") %>%
  mutate(ref = sub("_mobtyper", "", ref),
         ref = sub(".fasta_report.txt", "", ref)) %>%
  select(-file_id)

prokka_reports <- get_data(".",
                           pattern = "prokka_report",
                           convert = TRUE) %>%
  mutate(ref = sub("prokka_report_", "", ref),
         ref = sub(".txt", "", ref)) %>%
  rename("val" = `organism: Genus species strain `) %>%
  separate(val, c("key", "value"), ":") %>%
  spread(key, value) %>%
  select(ref, CDS, gene)

plasmidfinder_reports <- get_data(".",
                                  pattern = "plasfinder",
                                  convert = TRUE) %>%
  mutate(ref = sub("_plasfinder_results_tab.tsv", "", ref)) %>%
  group_by(ref) %>%
  summarise_all(list(func_paste)) %>%
  select(ref, Plasmid, Identity, Contig) %>%
  rename("pf_replicon" = Plasmid,
         "pf_replicon_identity" = Identity,
         "pf_contig" = Contig)

resfinder_reports <- get_data(".",
                              pattern = "resfinder",
                              convert = TRUE) %>%
  mutate(ref = sub("_resfinder_results_tab.tsv", "", ref)) %>%
  group_by(ref) %>%
  summarise_all(list(func_paste)) %>%
  select(ref, `Resistance gene`, Identity, Contig) %>%
  rename("rf_gene" = `Resistance gene`,
         "rf_gene_identity" = Identity,
         "rf_contig" = Contig)

virfinder_reports <- get_data(".",
                              pattern = "virfinder",
                              convert = TRUE) %>%
  mutate(ref = sub("_virfinder_results_tab.tsv", "", ref)) %>%
  group_by(ref) %>%
  summarise_all(list(func_paste)) %>%
  select(ref, `Virulence factor`, Identity, Contig) %>%
  rename("vf_gene" = `Virulence factor`,
         "vf_gene_identity" = Identity,
         "vf_contig" = Contig)

total_data <- prokka_reports %>%
  left_join(mobtyper_reports, by = "ref") %>%
  left_join(plasmidfinder_reports, by = "ref") %>%
  left_join(resfinder_reports, by = "ref") %>%
  left_join(virfinder_reports, by = "ref") %>%
  mutate(sample = sub("(.*?)_plasmid.+", "\\1", ref),
         plasmid = sub(".*?_(plasmid_.*?)", "\\1", ref)) %>%
  select(sample, plasmid, everything(), -ref)

write.table(total_data,
            "total_report.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

