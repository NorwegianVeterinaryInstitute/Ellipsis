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

