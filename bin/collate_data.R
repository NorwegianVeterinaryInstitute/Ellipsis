#!/usr/bin/env Rscript
# ABSTRACT: Collate data from Ellipsis into informative datasets
library(dplyr)
library(impoRt)
library(funtools)
library(tidyr)
library(readr)

path <- "."
args <- commandArgs(trailingOnly = TRUE)
run_ariba <- args[1]

# Functions
filter_flags <- function(df, acquired = TRUE) {
  acq_flags <- c(
    "27",
    "155",
    "411",
    "923",
    "795",
    "539",
    "667",
    "795"
  )
  
  int_flags <- c(
    "19",
    "27",
    "147",
    "155",
    "403",
    "411",
    "915",
    "923",
    "787",
    "795",
    "531",
    "539",
    "659",
    "667",
    "787",
    "795"
  )
  df <- df %>%
    select(ref, gene_names, flag, ref_ctg_change) %>%
    {if (acquired == TRUE) 
        filter(., flag %in% acq_flags) 
      else 
        filter(., flag %in% int_flags)} %>%
    select(-flag)
  
  stopifnot(dim(df)[1] != 0)
  
  return(df)
}

create_table <- function(df) {
  df %>%
    pivot_wider(names_from = gene_names, 
                values_from = ref_ctg_change,
                values_fn = func_paste) %>%
    pivot_longer(names_to = "gene",
                 values_to = "result",
                 cols = -ref) %>%
    mutate(
      result_total = ifelse(is.na(result), 0, 1),
      result_total = as.character(result_total)
    ) %>%
    select(-result) %>%
    rename("result" = result_total) %>%
    pivot_wider(names_from = gene,
                values_from = result,
                values_fn = func_paste)
}

fix_gene_names <- function(df, ending, db) {
  genes <- unique(df$ref_name)
  
  # correct the gene names with regex patterns
  # specific for each database
  if (db == "res") {
    new_names <- gsub("^(.*?)\\..*", "\\1", genes)
    new_names <- gsub("_", "", new_names, fixed = T)
    new_names <- gsub("-", "", new_names, fixed = T)
  }
  if (db == "vfdb") {
    new_names <- sub("\\..+", "", genes)
  }
  if (db == "virfinder") {
    new_names <- gsub("^(.+_[0-9]+)_.+", "\\1", genes)
  }
  
  # match the database names to the new names
  gene_names <- c()
  for (i in new_names) {
    p <- paste(tolower(substring(i, 1, 3)),
               substring(i, 4),
               sep = "",
               collapse = " ")
    gene_names <- c(gene_names, p)
  }
  df2 <- data.frame(genes, gene_names) %>%
    mutate(genes = as.character(genes)) %>%
    rename(ref_name = genes)
  
  df %>%
    left_join(df2, by = "ref_name") %>%
    mutate(
      gene_names = as.character(gene_names),
      ref = gsub(paste0("(.*?)", ending), "\\1", ref)
    )
}

# import data
contig_reports <- get_data(path,
                           pattern = "contig_report.txt",
                           convert = TRUE) %>%
  mutate(element = if_else(
           molecule_type == "chromosome",
           "chromosome",
           paste0(molecule_type, "_", primary_cluster_id)
           )) %>%
  select(sample_id, element, everything(), -c(ref, molecule_type)) %>%
  rename("ref" = sample_id)

mobtyper_reports <- get_data(path,
                             pattern = "mobtyper_results",
                             convert = TRUE) %>%
  mutate(element = paste0("plasmid_", sub(".+:(.*?)", "\\1", sample_id)),
         ref = sub("_mobtyper_results.txt", "", ref)) %>%
  select(ref, element, everything(), -sample_id)

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
  select(ref, element, CDS, gene, bases) %>%
  mutate(bases = gsub(" ", "", bases))

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
                              pattern = "resfinder_results_tab",
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
                              pattern = "virfinder_results_tab",
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

mlst_report <- get_data(path,
                        pattern = "assembly_mlst_report",
                        convert = TRUE) %>%
  mutate(FILE = sub(".fasta", "", FILE)) %>%
  select(-ref)

if (run_ariba == "true") {
  mlst_ariba <- get_data(path,
                         pattern = "ariba_mlst_report",
                         convert = TRUE) %>%
    mutate(ref = sub("_ariba_mlst_report.tsv","", ref))
  
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
  
  if (is.data.frame(df_vir) == TRUE & is.data.frame(df_res) == TRUE) {
    create_table(df_vir) %>%
      write_delim(path = "ariba_virfinder_results.txt", delim = "\t")
    create_table(df_res) %>%
      write_delim(file = "ariba_resfinder_results.txt", delim = "\t")
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
      write_delim(file = "ariba_resfinder_results.txt", delim = "\t")
  }
  
  if (is.data.frame(df_vir) == FALSE & is.data.frame(df_res) == FALSE) {
    "No genes passed quality checks" %>% 
      write_lines("ariba_virfinder_results.txt")
    "No genes passed quality checks" %>% 
      write_lines("ariba_resfinder_results.txt")
  }
} 

summary_report <- contig_reports %>%
  select(ref, element, size, gc, circularity_status) %>%
  left_join(prokka_reports[, c("ref", "element", "gene")], by = c("ref", "element")) %>%
  mutate(
    plasmid_size = ifelse(element != "chromosome", size, NA),
    plasmid_gene = ifelse(element != "chromosome", gene, NA)
  ) %>%
  mutate(plasmid_size = as.numeric(plasmid_size),
         plasmid_gene = as.numeric(plasmid_gene)) %>%
  group_by(ref) %>%
  mutate(plasmid = ifelse(element == "chromosome", 0, 1),
         n_plasmids = sum(plasmid)) %>%
  select(ref, circularity_status, n_plasmids) %>%
  summarise_all(list(func_paste)) %>%
  ungroup() %>%
  mutate(closed_genome = ifelse(grepl("incomplete", circularity_status), FALSE, TRUE)) %>%
  select(ref, closed_genome, everything(),-circularity_status) %>%
  left_join(mlst_report[,c("FILE","ST")], by = c("ref" = "FILE")) %>%
  rename("assembly_ST" = ST)

if (run_ariba == "true") {
  summary_report <- summary_report %>%
    left_join(mlst_ariba[,c("ref","ST")], by = "ref") %>%
    rename("reads_ST" = ST)
}

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

res_contig_truncated <- resfinder_reports %>%
  select(ref, element, rf_contig, rf_gene) %>%
  group_by(ref, element, rf_contig) %>%
  summarise_all(list(func_paste))

vir_contig_truncated <- virfinder_reports %>%
  select(ref, element, vf_contig, vf_gene) %>%
  group_by(ref, element, vf_contig) %>%
  summarise_all(list(func_paste))

plas_contig_truncated <- plasmidfinder_reports %>%
  select(ref, element, pf_contig, plasmidfinder_replicon) %>%
  group_by(ref, element, pf_contig) %>%
  summarise_all(list(func_paste))

contig_report <- contig_reports[, c(
  "ref",
  "element",
  "circularity_status",
  "contig_id",
  "size",
  "gc",
  "rep_type(s)",
  "relaxase_type(s)",
  "repetitive_dna_type",
  "filtering_reason"
)] %>%
  left_join(plas_contig_truncated, by = c("ref","element","contig_id" = "pf_contig")) %>%
  left_join(res_contig_truncated, by = c("ref","element","contig_id" = "rf_contig")) %>%
  left_join(vir_contig_truncated, by = c("ref","element","contig_id" = "vf_contig")) %>%
  select(-contig_id)

element_report <- contig_reports[, c("ref",
                                     "element",
                                     "circularity_status",
                                     "size",
                                     "rep_type(s)",
                                     "relaxase_type(s)")] %>%
  mutate(n = 1,
         size = as.numeric(size)) %>%
  group_by(ref, element) %>%
  mutate(total_size = sum(size),
         contigs = sum(n)) %>%
  na_if("-") %>%
  summarise_all(list(func_paste)) %>%
  mutate(closed = if_else(grepl("incomplete", circularity_status), FALSE, TRUE)) %>%
  ungroup() %>%
  select(
    ref,
    element,
    contigs,
    total_size,
    closed,
    everything(),
    -c(circularity_status, size, n)
  ) %>%
  left_join(plas_truncated, by = c("ref","element")) %>%
  left_join(res_truncated, by = c("ref", "element")) %>%
  left_join(vir_truncated, by = c("ref", "element")) %>%
  left_join(prokka_reports, by = c("ref","element")) %>%
  select(-bases) %>%
  relocate(CDS, .after = closed) %>%
  relocate(gene, .after = CDS)

# Write reports
write.table(summary_report,
            "summary_report.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

write.table(contig_report,
            "contig_report.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

write.table(element_report,
            "element_report.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

write.table(mobtyper_reports,
            "mobtyper_report.txt",
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



