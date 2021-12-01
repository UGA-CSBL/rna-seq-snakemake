# Snakemake pipeline for RNA-Seq analysis

## Overview

This pipeline runs the following steps:

1. Download mouse reference genome, transcriptome and gene annotation files from [GENCODE](https://www.gencodegenes.org/mouse/).
2. Run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://multiqc.info/) on the FASTQ files listed in the [design matrix](./config/design_matrix.csv).
3. Use [Salmon](https://salmon.readthedocs.io/en/stable/salmon.html) to quantify transcript abundance for each sample.
4. Use [STAR](https://github.com/alexdobin/STAR) to map reads to the genome into a BAM file, and count number of reads per gene while mapping.
5. Merge Salmon TPMs, Salmon read counts, and STAR read counts into three tables.
6. Use [DESeq2](https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html) to compare transcript abundance between samples.

## Getting started

First install Snakemake following the official [installation instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). Next, modify the [configuration file](./config/config.yaml) to set the reference file versions and FASTQ file locations.

Now we can simply run the pipeline with:

```bash
snakemake -c <num-of-cores> --use-conda
```

It's recommended to run the pipeline in `tmux` or `screen` as it can take a long time. Output to `stdout` are also saved to `.snakemake/log/`. Output of R scripts are saved to `logs/`.

This pipeline will genrate the following subdirectories in the `results/` directory:

-   `annotation/`: GENCODE reference files, Salmon index files, and STAR index files.
-   `qc/`: directory for FastQC and MultiQC output.
-   `salmon/`: directory for Salmon output.
-   `star/`: directory for STAR output.
-   `deseq2/`: directory for differential expression analysis results.
