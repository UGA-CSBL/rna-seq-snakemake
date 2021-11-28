library(tidyverse)

# Prepare file paths for the read count and TPM files -----
work_dir <- "~/data0/archive/JLU_mouse/2021_analysis"

design_mat <- data.table::fread(str_glue("{work_dir}/config/design_matrix.csv"))

s_ids <- unique(design_mat$Sample)
salmon_files <- str_glue("{work_dir}/quant/{s_ids}/quant.sf")
star_files <- str_glue("{work_dir}/star/{s_ids}/{s_ids}_ReadsPerGene.out.tab")

# Load transcript -> gene ID mapping from the GTF annotation file -----
tx2gene <- rtracklayer::import(
  str_glue("{work_dir}/annotation/gencode/gene_annotations.gencode"),
  format = "gtf"
)  # Or use GenomicFeatures::makeTxDbFromGFF
tx2gene <- tx2gene %>%
  as.data.frame() %>%
  dplyr::select(TXNAME = transcript_id, GENEID = gene_id) %>%
  distinct() %>%
  drop_na()

# Summarize Salmon transcript expression values at the gene level -----
salmon <- tximport::tximport(salmon_files, type = "salmon", tx2gene = tx2gene)

salmon_counts <- round(salmon$counts) %>%
  as.data.frame() %>%
  rename_with(~ s_ids) %>%
  rownames_to_column("Ensembl")

tpms <- salmon$abundance %>%
  as.data.frame() %>%
  rename_with(~ s_ids) %>%
  rownames_to_column("Ensembl")

# Concatenate STAR read counts into one table -----
star_counts <- data.frame(Ensembl = unique(tx2gene$GENEID))
for (i in seq_along(star_files)) {
  s_id <- s_ids[i]
  star <- data.table::fread(star_files[i]) %>%
    dplyr::select(Ensembl = V1, {{s_id}} := V2)

  star_counts <- full_join(star_counts, star, by = "Ensembl")
}

# Save files to disk -----
data.table::fwrite(salmon_counts, str_glue("{work_dir}/results/salmon_counts.csv"))
data.table::fwrite(tpms, str_glue("{work_dir}/results/TPMs.csv"))
data.table::fwrite(star_counts, str_glue("{work_dir}/results/star_counts.csv"))

# sessionInfo()
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.3 LTS
#
# Matrix products: default
# BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so
#
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
# [1] forcats_0.5.1   stringr_1.4.0   dplyr_1.0.7     purrr_0.3.4
# [5] readr_2.0.2     tidyr_1.1.4     tibble_3.1.5    ggplot2_3.3.5
# [9] tidyverse_1.3.1
#
# loaded via a namespace (and not attached):
# [1] MatrixGenerics_1.4.3        Biobase_2.52.0
# [3] httr_1.4.2                  bit64_4.0.5
# [5] vroom_1.5.5                 jsonlite_1.7.2
# [7] modelr_0.1.8                assertthat_0.2.1
# [9] stats4_4.1.1                GenomeInfoDbData_1.2.6
# [11] cellranger_1.1.0            Rsamtools_2.8.0
# [13] yaml_2.2.1                  pillar_1.6.4
# [15] backports_1.2.1             lattice_0.20-45
# [17] glue_1.4.2                  GenomicRanges_1.44.0
# [19] XVector_0.32.0              rvest_1.0.1
# [21] colorspace_2.0-2            Matrix_1.3-4
# [23] XML_3.99-0.8                pkgconfig_2.0.3
# [25] broom_0.7.9                 haven_2.4.3
# [27] zlibbioc_1.38.0             scales_1.1.1
# [29] tzdb_0.1.2                  BiocParallel_1.26.2
# [31] generics_0.1.1              IRanges_2.26.0
# [33] ellipsis_0.3.2              withr_2.4.2
# [35] SummarizedExperiment_1.22.0 BiocGenerics_0.38.0
# [37] cli_3.0.1                   magrittr_2.0.1
# [39] crayon_1.4.1                readxl_1.3.1
# [41] fs_1.5.0                    fansi_0.5.0
# [43] xml2_1.3.2                  tools_4.1.1
# [45] data.table_1.14.2           hms_1.1.1
# [47] BiocIO_1.2.0                lifecycle_1.0.1
# [49] matrixStats_0.61.0          S4Vectors_0.30.1
# [51] munsell_0.5.0               reprex_2.0.1
# [53] DelayedArray_0.18.0         Biostrings_2.60.2
# [55] compiler_4.1.1              GenomeInfoDb_1.28.4
# [57] rlang_0.4.12                grid_4.1.1
# [59] RCurl_1.98-1.5              tximport_1.20.0
# [61] rstudioapi_0.13             rjson_0.2.20
# [63] bitops_1.0-7                restfulr_0.0.13
# [65] gtable_0.3.0                DBI_1.1.1
# [67] R6_2.5.1                    GenomicAlignments_1.28.0
# [69] lubridate_1.7.10            rtracklayer_1.52.1
# [71] bit_4.0.4                   utf8_1.2.2
# [73] stringi_1.7.5               parallel_4.1.1
# [75] Rcpp_1.0.7                  vctrs_0.3.8
# [77] dbplyr_2.1.1                tidyselect_1.1.1
