log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(tidyverse)

# Prepare file paths for the read count and TPM files -----
star_files <- snakemake@input[["star_files"]]
tx2gene <- data.table::fread(snakemake@input[["tx2gene"]])

design_mat <- data.table::fread(snakemake@config[["samples"]])
s_ids <- unique(design_mat$Sample)

# Concatenate STAR read counts into one table -----
star_counts <- data.frame(Ensembl = unique(tx2gene$GENEID))
for (i in seq_along(star_files)) {
  s_id <- s_ids[i]
  star <- data.table::fread(star_files[i]) %>%
    dplyr::select(Ensembl = V1, {{s_id}} := V2)

  star_counts <- full_join(star_counts, star, by = "Ensembl")
}

# Move rows N_unmapped, N_ambiguous, N_noFeature, N_multimapping to the top -----
star_counts %>%
  arrange(desc(Ensembl)) %>%
  data.table::fwrite(snakemake@output[[1]])
