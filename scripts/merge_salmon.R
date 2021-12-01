log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(tidyverse)

# Summarize Salmon transcript expression values at the gene level -----
salmon_files <- snakemake@input[["salmon_files"]]
tx2gene <- data.table::fread(snakemake@input[["tx2gene"]])

# Extract column names from design matrix
design_mat <- data.table::fread(snakemake@config[["samples"]])
s_ids <- unique(design_mat$Sample)

salmon <- tximport::tximport(salmon_files, type = "salmon", tx2gene = tx2gene)

salmon_counts <- round(salmon$counts) %>%
  as.data.frame() %>%
  rename_with(~ s_ids) %>%
  rownames_to_column("Ensembl")

tpms <- salmon$abundance %>%
  as.data.frame() %>%
  rename_with(~ s_ids) %>%
  rownames_to_column("Ensembl")

# Save files to disk -----
data.table::fwrite(salmon_counts, snakemake@output[["counts"]])
data.table::fwrite(tpms, snakemake@output[["tpm"]])
