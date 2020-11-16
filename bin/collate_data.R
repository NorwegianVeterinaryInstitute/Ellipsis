#!/usr/bin/env Rscript
# ABSTRACT: Collate data from the plasmid pipeline into one table
library(dplyr)
library(impoRt)
library(funtools)
library(tidyr)
library(readr)

path <- "."
args <- commandArgs(trailingOnly = TRUE)
run_ariba <- args[1]

source("bin/collate_functions.R")

# import data
mobtyper_reports <- get_data(path,
                             pattern = "mobtyper",
                             convert = TRUE) %>%
  mutate(ref = sub("_mobtyper", "", ref),
         ref = sub(".fasta_report.txt", "", ref),
         element = sub(".*_(plasmid_.*)", "\\1", ref)) %>%
  rowwise() %>%
  mutate(ref = sub(paste0("_", element), "", ref)) %>%
  ungroup() %>%
  select(ref, element, everything(), -file_id)

prokka_reports <- get_data(path,
                           pattern = "prokka_report",
                           convert = TRUE) %>%
  mutate(ref = sub("prokka_report_", "", ref),
         ref = sub(".txt", "", ref),
         element = ifelse(grepl("chromosome", ref) == FALSE,
                          sub(".*_(plasmid_.'?)|", "\\1", ref), 
                          "chromosome")) %>%
  rowwise() %>%
  mutate(ref = sub(paste0("_", element), "", ref)) %>%
  ungroup() %>%
  rename("val" = `organism: Genus species strain `) %>%
  separate(val, c("key", "value"), ":") %>%
  spread(key, value) %>%
  select(ref, element, CDS, gene)

plasmidfinder_reports <- get_data(path,
                                  pattern = "plasfinder",
                                  convert = TRUE) %>%
  mutate(ref = sub("_plasfinder_results_tab.tsv", "", ref),
         Contig = sub(".*\\|", "", Contig),
         element = ifelse(grepl("chromosome", ref) == FALSE,
                          sub(".*_(plasmid_.'?)|", "\\1", ref), 
                          "chromosome")) %>%
  rowwise() %>%
  mutate(ref = sub(paste0("_", element), "", ref)) %>%
  ungroup() %>%
  select(ref, element, Plasmid, Identity, Contig) %>%
  rename("plasmidfinder_replicon" = Plasmid,
         "pf_replicon_identity" = Identity,
         "pf_contig" = Contig)

resfinder_reports <- get_data(path,
                              pattern = "resfinder",
                              convert = TRUE) %>%
  mutate(ref = sub("_resfinder_results_tab.tsv", "", ref),
         Contig = sub(".*\\|", "", Contig),
         element = ifelse(grepl("chromosome", ref) == FALSE,
                          sub(".*_(plasmid_.'?)|", "\\1", ref), 
                          "chromosome")) %>%
  rowwise() %>%
  mutate(ref = sub(paste0("_", element), "", ref)) %>%
  ungroup() %>%
  select(ref, element, `Resistance gene`, Identity, Contig) %>%
  rename("rf_gene" = `Resistance gene`,
         "rf_gene_identity" = Identity,
         "rf_contig" = Contig)

virfinder_reports <- get_data(path,
                              pattern = "virfinder",
                              convert = TRUE) %>%
  mutate(ref = sub("_virfinder_results_tab.tsv", "", ref),
         Contig = sub(".*\\|", "", Contig),
         element = ifelse(grepl("chromosome", ref) == FALSE,
                          sub(".*_(plasmid_.'?)|", "\\1", ref), 
                          "chromosome")) %>%
  rowwise() %>%
  mutate(ref = sub(paste0("_", element), "", ref)) %>%
  ungroup() %>%
  select(ref, element, `Virulence factor`, Identity, Contig) %>%
  rename("vf_gene" = `Virulence factor`,
         "vf_gene_identity" = Identity,
         "vf_contig" = Contig)

if (run_ariba == "true") {
<<<<<<< HEAD
  rawdata_res <- get_data(path,
                          pattern = "ariba_resfinder",
                          convert = TRUE) %>%
    fix_gene_names("_ariba_resfinder_report.tsv", db = "res")
  
  rawdata_vir <- get_data(path,
                          pattern = "ariba_virulence",
                          convert = TRUE) %>%
    fix_gene_names("_ariba_virulence_report.tsv", db = "virfinder")

  df_vir <- tryCatch(filter_flags(rawdata_vir), error = function(e) FALSE)
  df_res <- tryCatch(filter_flags(rawdata_res), error = function(e) FALSE)
=======
  ariba_res_report <- get_data(".",
                               pattern = "ariba_resfinder",
                               convert = TRUE) %>%
    fix_gene_names("_ariba_resfinder_report.tsv", db = "res") %>%
    create_table() %>%
    create_report() %>%
    rename("sample" = ref)
  
  ariba_vir_report <- get_data(".",
                               pattern = "ariba_virulence",
                               convert = TRUE) %>%
    fix_gene_names("_ariba_virulence_report.tsv", db = "virfinder") %>%
    create_table(acquired = FALSE) %>%
    create_report() %>%
    rename("sample" = ref)
>>>>>>> 82ebe5ac3a77802c4ac4543036aef818b36fc573
  
  if (is.data.frame(df_vir) == TRUE & is.data.frame(df_res) == TRUE) {
    create_table(df_vir) %>%
      write_delim(path = "ariba_virfinder_results.txt", delim = "\t")
    create_table(df_res) %>%
      write_delim(path = "ariba_resfinder_results.txt", delim = "\t")
  }
  
  if (is.data.frame(df_vir) == TRUE & is.data.frame(df_res) == FALSE) {
    create_table(df_vir) %>%
      write_delim(path = "ariba_virfinder_results.txt", delim = "\t")
    "No genes passed quality checks" %>% 
      write_lines("ariba_resfinder_results.txt")
  }
  
  if (is.data.frame(df_vir) == FALSE & is.data.frame(df_res) == TRUE) {
    "No genes passed quality checks" %>% 
      write_lines("ariba_virfinder_results.txt")
    create_table(df_res) %>%
      write_delim(path = "ariba_resfinder_results.txt", delim = "\t")
  }
  
  if (is.data.frame(df_vir) == FALSE & is.data.frame(df_res) == FALSE) {
    "No genes passed quality checks" %>% 
      write_lines("ariba_virfinder_results.txt")
    "No genes passed quality checks" %>% 
      write_lines("ariba_resfinder_results.txt")
  }
} 

summary_report <- mobtyper_reports %>%
    left_join(prokka_reports, by = c("ref", "element")) %>%
    mutate(
      total_length = as.numeric(total_length),
      num_contigs = as.numeric(num_contigs),
      gene = as.numeric(gene)
    ) %>%
    group_by(ref) %>%
    mutate(
      n_plasmids = length(element),
      size_range = paste(min(total_length), " - ", max(total_length)),
      contig_range = paste(min(num_contigs), " - ", max(num_contigs)),
      gene_range = paste(min(gene), " - ", max(gene))
    ) %>%
    select(ref, n_plasmids, size_range, contig_range, gene_range) %>%
    summarise_all(list(func_paste)) %>%
    rename("sample" = ref)  
  
# summarise and merge data
res_truncated <- resfinder_reports %>%
  select(ref, element, rf_gene) %>%
  group_by(ref, element) %>%
  summarise_all(list(func_paste))

vir_truncated <- virfinder_reports %>%
  select(ref, element, vf_gene) %>%
  group_by(ref, element) %>%
  summarise_all(list(func_paste))

plas_truncated <- plasmidfinder_reports %>%
  select(ref, element, plasmidfinder_replicon) %>%
  group_by(ref, element) %>%
  summarise_all(list(func_paste))

total_report <- prokka_reports %>%
  left_join(mobtyper_reports[, c(
    "ref",
    "element",
    "num_contigs",
    "total_length",
    "gc",
    "rep_type(s)",
    "PredictedMobility"
  )], by = c("ref", "element")) %>%
  left_join(plas_truncated, by = c("ref", "element")) %>%
  left_join(res_truncated, by = c("ref", "element")) %>%
  left_join(vir_truncated, by = c("ref", "element"))

# Write reports
write.table(summary_report,
            "summary_report.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

write.table(total_report,
            "total_report.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

write.table(resfinder_reports,
            "resfinder_report.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

write.table(virfinder_reports,
            "virfinder_report.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

write.table(plasmidfinder_reports,
            "plasmidfinder_report.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

