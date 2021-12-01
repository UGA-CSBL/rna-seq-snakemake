log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(DESeq2)
library(tidyverse)

# Configure multithread settings -----
parallel <- FALSE
if (snakemake@threads > 1) {
  library(BiocParallel)
  register(MulticoreParam(snakemake@threads))
  parallel <- TRUE
}

# Build design matrix for each timepoint -----
timepoint <- snakemake@wildcards[["tp"]]
software <- snakemake@wildcards[["software"]]

coldata <- data.table::fread(snakemake@config[["samples"]]) %>%
  mutate(Timepoint = str_sub(Sample, 2, 2)) %>%
  filter(Timepoint == timepoint) %>%
  select(SampleType, Sample) %>%
  arrange(SampleType) %>%
  filter(SampleType %in% c("Control", "Model")) %>%
  mutate(SampleType = as.factor(SampleType)) %>%
  column_to_rownames("Sample")

# Load STAR or Salmon read counts -----
cts <- data.table::fread(snakemake@input[[software]]) %>%
  select(Ensembl, all_of(rownames(coldata))) %>%
  filter(str_starts(Ensembl, "ENSMUSG")) %>%  # Drop N_* rows in STAR counts
  column_to_rownames("Ensembl")

# Perform differential expression analysis on each timepoint -----
message(str_glue("Running DESeq2 on {software} data for timepoint {timepoint}"))
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~ SampleType
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds$SampleType <- relevel(dds$SampleType, ref = "Control")

dds <- DESeq(dds, parallel = parallel)

results(dds, contrast = c("SampleType", "Model", "Control")) %>%
  as.data.frame() %>%
  rownames_to_column("Ensembl") %>%
  data.table::fwrite(snakemake@output[[1]])
