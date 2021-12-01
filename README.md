# Snakemake pipeline for RNA-Seq analysis

This pipeline runs the following steps:

1. Download mouse reference genome, transcriptome and gene annotation files from [GENCODE](https://www.gencodegenes.org/mouse/).
2. Run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://multiqc.info/) on the FASTQ files listed in the [design matrix](./config/design_matrix.csv).
3. Use [Salmon](https://salmon.readthedocs.io/en/stable/salmon.html) to quantify transcript abundance for each sample.
4. Use [STAR](https://github.com/alexdobin/STAR) to map reads to the genome into a BAM file, and count number of reads per gene while mapping. 
5. Merge Salmon TPMs, Salmon read counts, and STAR read counts into three tables.
6. Use [DESeq2](https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html) to compare transcript abundance between samples.

To run this pipeline, install and activate the conda environment with:

```bash
conda env create -f rnaseq.yml -n rnaseq
conda activate rnaseq
```

Then after modifying the [configuration](./config/config.yml) file, run the pipeline with:

```bash
snakemake -c <num-of-cores> --use-conda
```

This will genrate the following directories:

- `annotation/`: GENCODE reference files, Salmon index files, and STAR index files.
- `qc/`: directory for FastQC and MultiQC output.
- `quant/`: directory for Salmon output.
- `star/`: directory for STAR output.
