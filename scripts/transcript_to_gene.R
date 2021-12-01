log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(tidyverse)

annot_file <- snakemake@input[[1]]

# Alternatively, use GenomicFeatures::makeTxDbFromGFF
tx2gene <- rtracklayer::import(annot_file, format = "gtf")
tx2gene <- tx2gene %>%
  as.data.frame() %>%
  dplyr::select(TXNAME = transcript_id, GENEID = gene_id) %>%
  distinct() %>%
  drop_na()

data.table::fwrite(tx2gene, snakemake@output[[1]])
