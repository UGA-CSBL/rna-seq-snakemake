log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(biomaRt)

ensembl_version <- snakemake@params[["ensembl_version"]]
file_with_ensembl_ids <- snakemake@input[["ensembl_id_file"]]

ensembl <- useEnsembl(biomart = "genes", version = ensembl_version)
mouse <- useDataset("mmusculus_gene_ensembl", ensembl)
human <- useDataset("hsapiens_gene_ensembl", ensembl)

mouse_ids <- data.table::fread(file_with_ensembl_ids)$GENEID

id_mapping <- getLDS(
  mart = mouse,
  attributes = c("ensembl_gene_id_version", "external_gene_name"),
  martL = human,
  attributesL = c("ensembl_gene_id_version", "external_gene_name"),
  filters = "ensembl_gene_id_version",
  values = mouse_ids,
  uniqueRows = T
)

colnames(id_mapping) <- c("mouse_ensembl", "mouse_symbol",
                          "human_ensembl", "human_symbol")

data.table::fwrite(id_mapping, snakemake@output[[1]])
